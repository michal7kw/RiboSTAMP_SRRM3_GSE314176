#!/bin/bash
# =============================================================================
# pc_pipeline/04 — pbmm2 align segmented S-reads to mm10.
# Sequential per SRR; each uses NUM_THREADS.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

if [ ! -f "${MM10_REF}" ]; then
  echo "ERROR: mm10 reference not found at ${MM10_REF}." >&2
  echo "Run pc_pipeline/setup_references.sh first." >&2
  exit 1
fi

for SRR in ${SRR_LIST}; do
  IN_BAM="${SPLIT_DIR}/${SRR}.segmented.bam"
  OUT_BAM="${ALIGN_DIR}/${SRR}.mm10.bam"

  if [ ! -f "${IN_BAM}" ]; then
    echo "ERROR: missing ${IN_BAM}; run 03_skera_split.sh first." >&2
    exit 1
  fi
  if [ -s "${OUT_BAM}" ]; then
    log_step "${SRR}: already aligned; skipping."
    continue
  fi

  log_step "${SRR}: pbmm2 align (${NUM_THREADS} threads)"
  pbmm2 align \
      --preset ISOSEQ \
      --sort \
      --num-threads "${NUM_THREADS}" \
      "${MM10_REF}" \
      "${IN_BAM}" \
      "${OUT_BAM}"

  samtools index -@ 4 "${OUT_BAM}"

  total=$(samtools view -c "${OUT_BAM}")
  mapped=$(samtools view -c -F 4 "${OUT_BAM}")
  rate=$(awk -v m="${mapped}" -v t="${total}" 'BEGIN{printf "%.2f", 100*m/t}')
  echo "${SRR}: total=${total} mapped=${mapped} rate=${rate}%"

  if (( $(echo "${rate} < 80" | bc -l) )); then
    echo "WARNING: ${SRR} mapping rate ${rate}% < 80%. Investigate before continuing." >&2
  fi
done

log_step "All alignments done"
ls -lh "${ALIGN_DIR}"/*.mm10.bam
