#!/bin/bash
# =============================================================================
# pc_pipeline/06 — ON-DEMAND full-transcriptome IsoQuant.
# Submit only if 05 found a novel isoform and you want to confirm
# locus-subsetting didn't bias read support. ~24 h on this PC.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

log_step "Full-transcriptome IsoQuant (this is the long one — ~24 h)"
isoquant \
    --reference "${MM10_REF}" \
    --genedb "${GENCODE_GTF}" \
    --bam "${ALIGN_DIR}"/*.mm10.bam \
    --data_type pacbio_ccs \
    --model_construction_strategy sensitive_pacbio \
    --genedb_output \
    --output "${ISOQ_FULL_DIR}" \
    --threads "${NUM_THREADS}" \
    --sqanti_output \
    --count_exons \
    --read_group tag:CB

log_step "Done"
