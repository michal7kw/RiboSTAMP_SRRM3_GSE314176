#!/bin/bash
# =============================================================================
# 02 — two-phase alignment: fast minimap2 pre-filter, then pbmm2 splice-aware.
#
# Why two phases:
#   pbmm2 with ISOSEQ preset processes ~100 MB/min regardless of reference
#   size — that's ~10 hours per 60 GB BAM. The bottleneck is per-read
#   splice-aware DP. But Srrm3 reads are <0.1% of the input.
#
# Phase A: minimap2 with -x map-pb (NO splice) is ~10× faster as a coarse
#   filter. Most reads either seed-match the locus (= candidates) or get
#   rejected at the seed stage (rejected reads cost almost no CPU).
#
# Phase B: pbmm2 with ISOSEQ on the small candidate set is fast because
#   the input is now ~few MB instead of ~60 GB.
#
# Empirical expected runtime: ~5-10 min Phase A per BAM + ~1 min Phase B
#   = ~30 min total for all 3 SRRs (vs ~30 hours with the naive approach).
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

if [ ! -f "${LOCUS_FA}" ]; then
  echo "ERROR: ${LOCUS_FA} not found. Run 01_extract_locus.sh first." >&2
  exit 1
fi

CANDIDATE_DIR="${TGT_DIR}/candidate_reads"
mkdir -p "${CANDIDATE_DIR}"

SKERA_DIR="${BASE_DIR}/skera_split"
shopt -s nullglob
SEGMENTED_BAMS=("${SKERA_DIR}"/SRR*.segmented.bam)
shopt -u nullglob

if [ "${#SEGMENTED_BAMS[@]}" -eq 0 ]; then
  echo "ERROR: no segmented BAMs in ${SKERA_DIR}." >&2
  exit 1
fi

for IN_BAM in "${SEGMENTED_BAMS[@]}"; do
  NAME="$(basename "${IN_BAM}" .segmented.bam)"
  CANDIDATE_BAM="${CANDIDATE_DIR}/${NAME}.candidates.bam"
  CANDIDATE_FQ="${CANDIDATE_DIR}/${NAME}.candidates.fastq.gz"
  OUT_BAM="${TGT_ALIGN_DIR}/${NAME}.locus.bam"
  SRR_LOG="${TGT_LOG_DIR}/02_align_${NAME}.log"

  if [ -s "${OUT_BAM}" ]; then
    log_step "${NAME}: already aligned to locus; skipping."
    continue
  fi

  {
    # ---------------------------------------------------------------------
    # Phase A — minimap2 -x map-pb pre-filter (NO splice; fast)
    # ---------------------------------------------------------------------
    # NOTE: filter tightened after first run produced 3M candidate reads
    # per BAM (vs target ~10-100K real Srrm3 reads). Phase B then took
    # 4.5 hours per BAM instead of expected 5-15 min. The fix:
    #   --secondary=no    (no secondary alignments output)
    #   -p 0.5            (require >=50% of best DP score for primary —
    #                      rejects reads with marginal seed matches)
    #   -q 30 in samtools (require MAPQ >= 30 — drops mappings to repeats
    #                      and similar low-confidence alignments)
    # Combined with the 30kb (vs 200kb) locus FASTA, this should give
    # ~10-50× fewer candidates. See LEARNING_07.
    if [ ! -s "${CANDIDATE_FQ}" ]; then
      echo "[$(date '+%H:%M:%S')] ${NAME} — Phase A: minimap2 pre-filter (map-pb, strict)"
      samtools fastq -@ 4 "${IN_BAM}" 2>/dev/null \
        | minimap2 -t "${TGT_NUM_THREADS}" \
                   -ax map-pb \
                   --secondary=no \
                   -p 0.5 \
                   "${LOCUS_FA}" - 2>/dev/null \
        | samtools view -b -F 4 -q 30 -@ 2 - 2>/dev/null \
        | samtools fastq -@ 2 - 2>/dev/null \
        | gzip -1 > "${CANDIDATE_FQ}.tmp"
      mv "${CANDIDATE_FQ}.tmp" "${CANDIDATE_FQ}"
    fi

    n_candidates=$(zcat "${CANDIDATE_FQ}" | awk 'NR%4==1' | wc -l)
    echo "[$(date '+%H:%M:%S')] ${NAME}: ${n_candidates} candidate reads ($(stat -c%s "${CANDIDATE_FQ}") bytes gz)"

    if [ "${n_candidates}" -lt 1 ]; then
      echo "WARNING: ${NAME} has 0 candidate reads — writing empty output." >&2
      # Create empty BAM with valid header so downstream scripts don't crash
      samtools view -bH /dev/null > "${OUT_BAM}" 2>/dev/null || touch "${OUT_BAM}"
      continue
    fi

    # ---------------------------------------------------------------------
    # Phase B — pbmm2 splice-aware alignment of candidates
    # ---------------------------------------------------------------------
    echo "[$(date '+%H:%M:%S')] ${NAME} — Phase B: pbmm2 ISOSEQ on candidates"
    pbmm2 align \
        --preset ISOSEQ \
        --sort \
        --num-threads "${TGT_NUM_THREADS}" \
        "${LOCUS_FA}" \
        "${CANDIDATE_FQ}" \
      | samtools view -@ 2 -b -h -F 4 - > "${OUT_BAM}"

    samtools index -@ 2 "${OUT_BAM}"

    bytes=$(stat -c%s "${OUT_BAM}")
    reads=$(samtools view -c "${OUT_BAM}")
    echo "[$(date '+%H:%M:%S')] ${NAME}: ${reads} aligned reads (${bytes} bytes BAM)"
  } 2>&1 | tee "${SRR_LOG}"
done

log_step "Done"
ls -lh "${TGT_ALIGN_DIR}/"*.bam 2>/dev/null
echo
echo "=== candidate read counts (Phase A output) ==="
for f in "${CANDIDATE_DIR}"/*.candidates.fastq.gz; do
  [ -f "$f" ] || continue
  n=$(zcat "$f" | awk 'NR%4==1' | wc -l)
  printf "  %s: %d candidates\n" "$(basename $f)" "$n"
done
