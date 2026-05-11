#!/bin/bash
# =============================================================================
# 08_subset_locus.sh — subset full STARsolo BAMs to the Srrm3 locus
#
# The full STARsolo BAMs are ~50-80 GB each. For per-cell junction PSI we
# only need reads at the cassette locus (chr5:135,680,000-135,960,000 mm10),
# which is well under 1 GB per sample. Subsetting first means:
#   - junction calc reads tiny BAMs in seconds, not minutes
#   - we can keep small subsets long-term and free disk on the full BAMs
#
# Resume-friendly: skip if {GSM}_srrm3.bam already exists.
#
# Outputs:
#   ${SR_BAM_DIR_10X}/{GSM}_srrm3.bam      (locus subset, sorted, CB-tagged)
#   ${SR_BAM_DIR_10X}/{GSM}_srrm3.bam.bai
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/08_subset_locus_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "08_subset_locus.sh — subset full BAMs to Srrm3 locus"

activate_env
which samtools >/dev/null

LOCUS_BED="${THIS_DIR}/../../hpc/srrm3_locus.bed"
if [[ ! -s "${LOCUS_BED}" ]]; then
    echo "ERROR: missing ${LOCUS_BED}"
    exit 1
fi
echo "Locus BED: ${LOCUS_BED}"
cat "${LOCUS_BED}"

for GSM in "${SR_10X_GSMS[@]}"; do
    FULL="${SR_BAM_DIR_10X}/${GSM}.bam"
    SUBSET="${SR_BAM_DIR_10X}/${GSM}_srrm3.bam"

    if [[ ! -s "${FULL}" ]]; then
        echo "WARN: full BAM missing for ${GSM} (${FULL}); skip"
        continue
    fi
    if [[ -s "${SUBSET}" ]]; then
        echo "[skip] ${SUBSET} exists ($(du -h "${SUBSET}" | cut -f1))"
        continue
    fi

    echo ""
    echo "Subsetting ${GSM}: ${FULL} → ${SUBSET}"
    samtools view -bh -L "${LOCUS_BED}" -@ "${SR_NUM_THREADS}" "${FULL}" > "${SUBSET}"
    samtools index -@ "${SR_NUM_THREADS}" "${SUBSET}"
    echo "  ${GSM}: $(du -h "${SUBSET}" | cut -f1)  $(samtools view -c "${SUBSET}") reads"
done

log_step "08_subset_locus.sh complete"
echo ""
ls -lh "${SR_BAM_DIR_10X}"/*_srrm3.bam 2>&1
echo ""
echo "Next: bash 09_per_cell_progenitor_psi.sh"
