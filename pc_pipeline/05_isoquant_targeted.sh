#!/bin/bash
# =============================================================================
# pc_pipeline/05 — IsoQuant on Srrm3 locus subset, sensitivity raised for
# known-target novel-isoform discovery (see LEARNING_06).
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

if [ ! -f "${LOCUS_BED}" ]; then
  echo "ERROR: ${LOCUS_BED} not found." >&2
  exit 1
fi

log_step "Step 1 — subset all aligned BAMs to Srrm3 locus"
SUBSET_BAMS=()
for BAM in "${ALIGN_DIR}"/*.mm10.bam; do
  [ -f "${BAM}" ] || continue
  NAME="$(basename "${BAM}" .mm10.bam)"
  OUT_BAM="${SUBSET_DIR}/${NAME}.srrm3.bam"
  if [ ! -s "${OUT_BAM}" ]; then
    samtools view -@ "${NUM_THREADS}" -b -L "${LOCUS_BED}" "${BAM}" > "${OUT_BAM}"
    samtools index -@ 4 "${OUT_BAM}"
  fi
  SUBSET_BAMS+=("${OUT_BAM}")
done
echo "Subsetted ${#SUBSET_BAMS[@]} BAMs."

if [ "${#SUBSET_BAMS[@]}" -eq 0 ]; then
  echo "ERROR: no aligned BAMs in ${ALIGN_DIR}." >&2
  exit 1
fi

log_step "Step 2 — IsoQuant (sensitive_pacbio mode, novel discovery enabled)"
isoquant \
    --reference "${MM10_REF}" \
    --genedb "${GENCODE_GTF}" \
    --bam "${SUBSET_BAMS[@]}" \
    --data_type pacbio_ccs \
    --model_construction_strategy sensitive_pacbio \
    --output "${ISOQ_TGT_DIR}" \
    --threads "${NUM_THREADS}" \
    --sqanti_output \
    --count_exons

log_step "Done"
ls -lh "${ISOQ_TGT_DIR}"
