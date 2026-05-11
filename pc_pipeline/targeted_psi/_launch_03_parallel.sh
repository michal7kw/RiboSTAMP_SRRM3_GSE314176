#!/bin/bash
# Detach _run_03_parallel.sh as a background process.
set -e
HERE="$(dirname "$(readlink -f "$0")")"
source "${HERE}/00_config.sh"

LAUNCH_LOG="${TGT_LOG_DIR}/master_03_parallel_$(date +%Y%m%d_%H%M%S).log"
INFO_FILE="${TGT_LOG_DIR}/last_run.info"

source "${CONDA_BASE}/etc/profile.d/conda.sh"

nohup bash "${HERE}/_run_03_parallel.sh" > "${LAUNCH_LOG}" 2>&1 &
PID=$!
disown

{
  echo "PID=${PID}"
  echo "LOG=${LAUNCH_LOG}"
  echo "STARTED=$(date -Iseconds)"
} > "${INFO_FILE}"

echo "Step 03 (parallel) + step 04 launched."
echo "  PID: ${PID}"
echo "  Log: ${LAUNCH_LOG}"
