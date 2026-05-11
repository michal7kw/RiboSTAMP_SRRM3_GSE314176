#!/bin/bash
# =============================================================================
# pc_pipeline/04_fast_align_srrm3.sh
#
# Fast-path alternative to 04_align_pbmm2.sh for Path A (de-novo Srrm3 isoform
# discovery only). Reuses the minimap2-prefiltered candidate FASTQs already
# produced by targeted_psi/02_align_to_locus.sh — those are a small Srrm3-
# enriched subset (~2M reads/SRR vs ~50M reads in the full segmented BAMs),
# so realigning them to the full mm10 genome takes ~10-30 min/SRR instead of
# the ~6-8 h that naive whole-genome pbmm2 needs.
#
# Outputs land in ALIGN_DIR/<SRR>.mm10.bam — the same path the canonical
# 05_isoquant_targeted.sh expects, so that script runs unchanged afterwards.
#
# Trade-offs vs the canonical 04:
#   - LOSES the CB cell-barcode tag (the candidate FASTQs went through
#     `samtools fastq` which strips BAM tags). This is fine for STRUCTURAL
#     isoform discovery; if per-cluster expression of any novel isoform is
#     wanted later, we can recover CB by joining read IDs back to the
#     segmented BAM.
#   - MAY MISS Srrm3 reads whose only mm10 anchor is in the 5' half of
#     Srrm3 (chr5:135,806,890-135,839,721). Reads spanning the full
#     transcript still pass the prefilter via 3' seed matches, so this is
#     unlikely to lose much. Coverage is verified at the end of this script.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

CAND_DIR="${BASE_DIR}/targeted_psi/candidate_reads"

if [ ! -f "${MM10_REF}" ]; then
  echo "ERROR: mm10 reference not found at ${MM10_REF}." >&2
  exit 1
fi

if [ ! -d "${CAND_DIR}" ]; then
  echo "ERROR: candidate reads directory not found: ${CAND_DIR}" >&2
  echo "Run pc_pipeline/targeted_psi/02_align_to_locus.sh first." >&2
  exit 1
fi

shopt -s nullglob
CAND_FQS=("${CAND_DIR}"/*.candidates.fastq.gz)
shopt -u nullglob

if [ "${#CAND_FQS[@]}" -eq 0 ]; then
  echo "ERROR: no candidate FASTQs in ${CAND_DIR}." >&2
  exit 1
fi

log_step "Fast-path Srrm3 alignment — ${#CAND_FQS[@]} candidate FASTQs → mm10"

for FQ in "${CAND_FQS[@]}"; do
  NAME="$(basename "${FQ}" .candidates.fastq.gz)"
  OUT_BAM="${ALIGN_DIR}/${NAME}.mm10.bam"

  if [ -s "${OUT_BAM}" ] && [ -s "${OUT_BAM}.bai" ]; then
    log_step "${NAME}: already aligned to mm10; skipping."
    continue
  fi

  log_step "${NAME}: pbmm2 ISOSEQ → mm10 (${NUM_THREADS} threads)"
  pbmm2 align \
      --preset ISOSEQ \
      --sort \
      --num-threads "${NUM_THREADS}" \
      "${MM10_REF}" \
      "${FQ}" \
      "${OUT_BAM}"

  samtools index -@ 4 "${OUT_BAM}"

  total=$(samtools view -c "${OUT_BAM}")
  mapped=$(samtools view -c -F 4 "${OUT_BAM}")
  if [ "${total}" -gt 0 ]; then
    rate=$(awk -v m="${mapped}" -v t="${total}" 'BEGIN{printf "%.2f", 100*m/t}')
  else
    rate="0.00"
  fi
  echo "${NAME}: total=${total} mapped=${mapped} rate=${rate}%"
done

log_step "Coverage check — Srrm3 5' half vs 3' half"
# Srrm3 mm10: chr5:135,806,890-135,874,772 (gene span from Gencode vM25)
# 5' half: chr5:135,806,890-135,840,831 (~34 kb)
# 3' half: chr5:135,840,831-135,874,772 (~34 kb)
# If the prefilter missed 5' reads, the 5' count will be much lower than 3'.
for BAM in "${ALIGN_DIR}"/*.mm10.bam; do
  [ -f "${BAM}" ] || continue
  NAME="$(basename "${BAM}" .mm10.bam)"
  n5=$(samtools view -c "${BAM}" chr5:135806890-135840831)
  n3=$(samtools view -c "${BAM}" chr5:135840831-135874772)
  printf "  %s: 5'=%d  3'=%d\n" "${NAME}" "${n5}" "${n3}"
done

log_step "Done"
ls -lh "${ALIGN_DIR}"/*.mm10.bam
