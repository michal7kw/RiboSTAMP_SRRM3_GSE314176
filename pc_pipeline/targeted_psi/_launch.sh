#!/bin/bash
# Detach run_all.sh as a background process. Mirrors pc_pipeline/_launch_full.sh.
set -e

HERE="$(dirname "$(readlink -f "$0")")"
source "${HERE}/00_config.sh"

source "${CONDA_BASE}/etc/profile.d/conda.sh"

LAUNCH_LOG="${TGT_LOG_DIR}/master_$(date +%Y%m%d_%H%M%S).log"
INFO_FILE="${TGT_LOG_DIR}/last_run.info"

nohup bash "${HERE}/run_all.sh" > "${LAUNCH_LOG}" 2>&1 &
PID=$!
disown

{
  echo "PID=${PID}"
  echo "LOG=${LAUNCH_LOG}"
  echo "STARTED=$(date -Iseconds)"
} > "${INFO_FILE}"

echo "Targeted PSI pipeline launched."
echo "  PID:  ${PID}"
echo "  Log:  ${LAUNCH_LOG}"
echo "  Info: ${INFO_FILE}"
