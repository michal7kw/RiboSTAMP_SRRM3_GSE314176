#!/bin/bash
# =============================================================================
# One-time provisioning of the longread_iso conda env.
#
# Runs interactively on a login node — DO NOT submit via sbatch.
# Wall time ~10 min for the mamba solve + download + link step.
# Idempotent: if the env already exists, exits 0 without doing anything.
# =============================================================================

set -euo pipefail

ENV_PREFIX="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso"
ENV_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/RiboSTAMP_SRRM3_GSE314176/env/longread_iso.yml"

if [[ -d "${ENV_PREFIX}" ]]; then
    echo "Env already exists at ${ENV_PREFIX} — nothing to do."
    exit 0
fi

if [[ ! -f "${ENV_FILE}" ]]; then
    echo "ERROR: env file not found at ${ENV_FILE}" >&2
    exit 1
fi

echo "============================================"
echo "Creating longread_iso env"
echo "Prefix: ${ENV_PREFIX}"
echo "Spec:   ${ENV_FILE}"
echo "Start:  $(date)"
echo "============================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
mamba env create --prefix "${ENV_PREFIX}" --file "${ENV_FILE}"

conda activate "${ENV_PREFIX}"
python -c "import pysam, pandas, scipy; print('OK: pysam, pandas, scipy importable')"
which python samtools

echo ""
echo "=== Done: $(date) ==="
