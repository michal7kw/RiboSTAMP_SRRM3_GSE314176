#!/bin/bash
# =============================================================================
# pc_pipeline/02 — fetch PacBio HiFi BAMs ready for skera split.
#
# THE PROBLEM:
#   NCBI's SRA-Normalized format (what `prefetch` downloads) STRIPS PacBio
#   aux tags (qs/pw/ip/np/rq) that skera + IsoQuant need. The truly original
#   BAMs ARE preserved at AWS S3 (s3://sra-pub-src-7/...) but in
#   REQUESTER-PAYS buckets (~$16 USD egress for 3 BAMs), needing AWS creds.
#
# THE SOLUTION (this script):
#   1. prefetch the .sra (free, anonymous)
#   2. sam-dump to SAM (preserves SEQ + per-base QUAL)
#   3. SYNTHESIZE the few aux fields skera demands:
#        - PacBio-format QNAME: <movie>/<zmw>/ccs
#        - @RG header line + matching RG:Z: per-record tag
#        - rq:f:0.99 (predicted read accuracy — HiFi by definition is high-Q)
#        - np:i:15   (number of passes — typical for Revio HiFi)
#        - qs:i:0  qe:i:<seq_len>  (HQ region = entire read)
#   4. Pipe to samtools view -b -> compressed BAM
#
# Critically: kinetic tags (pw/ip) that CANNOT be synthesized are also NOT
# needed by skera, pbmm2, or IsoQuant in our pipeline. Verified empirically:
# 10K test reads -> 158K S-reads at 15.8 segments/read, 94.8% full arrays
# (exactly what MAS-16 expects).
#
# See docs/LEARNING_05_masseq_skera.md for why we don't use the AWS path.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

# Synth defaults — applied to every record. These have to look like a
# PacBio Revio HiFi run to satisfy pbbam's validator, but the values are
# nominal because skera doesn't filter on them.
SYNTH_RG_ID="m_synthetic_pacbio_RG"
SYNTH_RG_LINE="@RG	ID:${SYNTH_RG_ID}	PL:PACBIO	PM:REVIO	PU:${SYNTH_RG_ID}	SM:1	DS:READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100"

for SRR in ${SRR_LIST}; do
  OUT_BAM="${RAW_DIR}/${SRR}.unaligned.bam"
  if [ -s "${OUT_BAM}" ]; then
    log_step "${SRR}: BAM already exists; skipping. Delete it to re-process."
    continue
  fi

  # Step A — prefetch (idempotent; resumes partials)
  log_step "${SRR}: prefetch (~38 GB SRA-Normalized)"
  prefetch --output-directory "${RAW_DIR}" --max-size u --resume yes "${SRR}"

  SRA_PATH="${RAW_DIR}/${SRR}/${SRR}.sra"
  [ -f "${SRA_PATH}" ] || SRA_PATH="${RAW_DIR}/${SRR}.sra"
  if [ ! -f "${SRA_PATH}" ]; then
    echo "ERROR: ${SRR}.sra not where expected." >&2
    find "${RAW_DIR}" -name "${SRR}*.sra" >&2
    exit 1
  fi

  # Step B — sam-dump | synth-tags awk | samtools view -b
  log_step "${SRR}: sam-dump + synthesize PacBio aux fields -> BAM"
  {
    printf "@HD\tVN:1.6\tSO:unknown\tpb:5.0.0\n"
    printf "%s\n" "${SYNTH_RG_LINE}"
    sam-dump --unaligned "${SRA_PATH}" 2>/dev/null \
      | awk -v rg="${SYNTH_RG_ID}" 'BEGIN{OFS="\t"} {
          if ($0 ~ /^@/) next   # drop sam-dump headers; we already wrote our own
          seq_len = length($10)
          # Build a PacBio-style QNAME: <movie>/<zmw>/ccs
          $1 = rg "/" NR "/ccs"
          # Reconstruct the standard SAM fields 1..11 cleanly, dropping
          # whatever aux tags sam-dump emitted (only RG in practice).
          out = $1
          for (i = 2; i <= 11; i++) out = out "\t" $i
          # Append our standardized aux tags
          printf "%s\tRG:Z:%s\trq:f:0.99\tnp:i:15\tqs:i:0\tqe:i:%d\n", out, rg, seq_len
        }'
  } | samtools view -@ "${NUM_THREADS}" -b -h - > "${OUT_BAM}"

  bytes=$(stat -c%s "${OUT_BAM}")
  if [ "${bytes}" -lt 1000000000 ]; then
    echo "WARNING: ${OUT_BAM} is only ${bytes} bytes — likely incomplete." >&2
  fi
  log_step "${SRR}: $(numfmt --to=iec ${bytes}) BAM written."

  # Free the .sra cache once we have the BAM (recovers ~36 GB per run)
  rm -rfv "${RAW_DIR}/${SRR}/" 2>/dev/null || true
done

echo ""
log_step "All runs done"
ls -lh "${RAW_DIR}"/*.unaligned.bam

# Sanity check on the synth tags. We disable pipefail locally — samtools | head
# triggers SIGPIPE on samtools which would otherwise bubble up as a failure
# even though the data is fine.
set +eo pipefail
echo ""
echo "=== sanity check: aux tags + QNAME format ==="
for BAM in "${RAW_DIR}"/SRR*.unaligned.bam; do
  echo -n "$(basename ${BAM}): "
  first=$(samtools view "${BAM}" 2>/dev/null | head -1 2>/dev/null)
  qname=$(printf '%s' "${first}" | cut -f1)
  ntags=$(printf '%s' "${first}" | awk -F'\t' '{c=0; for(i=12;i<=NF;i++) c++; print c}')
  echo "QNAME=${qname:0:50}... ${ntags} aux tags"
done
set -eo pipefail
