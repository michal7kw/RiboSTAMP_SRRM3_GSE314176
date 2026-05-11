#!/bin/bash
# =============================================================================
# pc_pipeline/07 — SQANTI3 classification (FSM/ISM/NIC/NNC + artifact flags).
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
# SQANTI3 has its own conda env (not in bioconda). See setup_sqanti3.sh.
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SQANTI3_ENV}"

if [ ! -d "${SQANTI3_HOME}" ]; then
  echo "ERROR: SQANTI3 not found at ${SQANTI3_HOME}." >&2
  echo "Run pc_pipeline/setup_sqanti3.sh first." >&2
  exit 1
fi

NOVEL_GTF="${ISOQ_TGT_DIR}/OUT/OUT.transcript_models.gtf"
if [ ! -f "${NOVEL_GTF}" ]; then
  echo "ERROR: ${NOVEL_GTF} not found; run 05_isoquant_targeted.sh first." >&2
  exit 1
fi

log_step "SQANTI3 classification (using ${SQANTI3_HOME})"
# SQANTI3 v5.4 hardcodes some relative log paths under its own install dir,
# so invoke it from SQANTI3_HOME and use absolute paths for I/O.
( cd "${SQANTI3_HOME}" && \
  python "${SQANTI3_HOME}/sqanti3_qc.py" \
      --isoforms "${NOVEL_GTF}" \
      --refGTF "${GENCODE_GTF}" \
      --refFasta "${MM10_REF}" \
      -o ribostamp_srrm3 \
      -d "${SQANTI_DIR}" \
      --cpus "${NUM_THREADS}" \
      --report skip )

log_step "Done"
echo "Key outputs:"
echo "  ${SQANTI_DIR}/ribostamp_srrm3_classification.txt"
echo "  ${SQANTI_DIR}/ribostamp_srrm3_corrected.gtf"
