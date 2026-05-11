#!/bin/bash
# =============================================================================
# Re-quantify the same BAMs against the extended annotation (reference +
# discovered novel transcripts) so that OUT.read_assignments.tsv.gz contains
# per-read assignments to the NOVEL transcripts as well.
#
# This is a strictly assignment+quantification pass (no model construction),
# so it should run in seconds.
# =============================================================================
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

EXT_GTF="${ISOQ_TGT_DIR}/OUT/OUT.extended_annotation.gtf"
OUT_DIR="${ISOQ_TGT_DIR}_quant"
mkdir -p "${OUT_DIR}"

if [ ! -f "${EXT_GTF}" ]; then
  echo "ERROR: ${EXT_GTF} not found. Run 05_isoquant_targeted.sh first." >&2
  exit 1
fi

shopt -s nullglob
SUBSET_BAMS=("${SUBSET_DIR}"/*.srrm3.bam)
shopt -u nullglob
if [ "${#SUBSET_BAMS[@]}" -eq 0 ]; then
  echo "ERROR: no Srrm3 subset BAMs at ${SUBSET_DIR}" >&2
  exit 1
fi

log_step "IsoQuant re-quantify against extended annotation (no_model_construction)"
isoquant \
    --reference "${MM10_REF}" \
    --genedb "${EXT_GTF}" \
    --bam "${SUBSET_BAMS[@]}" \
    --data_type pacbio_ccs \
    --no_model_construction \
    --output "${OUT_DIR}" \
    --threads "${NUM_THREADS}" \
    --sqanti_output \
    --count_exons

log_step "Done"
ls -lh "${OUT_DIR}"
