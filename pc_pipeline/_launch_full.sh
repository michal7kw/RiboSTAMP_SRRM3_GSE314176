#!/bin/bash
# Internal helper — invoked once to detach the full pipeline as a background
# process. Don't call this directly; it's used by the WSL launch flow so we
# avoid quoting hell with nohup + & through `wsl bash -lc`.

set -e

# Resolve paths relative to this script's location
HERE="$(dirname "$(readlink -f "$0")")"
REPO_ROOT="$(dirname "${HERE}")"

# Sync repo
cd "${REPO_ROOT}"
git pull --ff-only 2>&1 | tail -3

# Make scripts executable in case git didn't preserve it on Windows side
chmod +x pc_pipeline/*.sh

# Ensure data + logs dirs exist
DATA_BASE="/mnt/e/RiboSTAMP_SRRM3_GSE314176/data"
mkdir -p "${DATA_BASE}/logs"

# Initialize conda for the detached subshell
source "${HOME}/miniconda3/etc/profile.d/conda.sh"

# Compute the log + info file paths
LAUNCH_LOG="${DATA_BASE}/logs/master_$(date +%Y%m%d_%H%M%S).log"
INFO_FILE="${DATA_BASE}/logs/last_run.info"

# Detach the actual pipeline run
nohup bash "${HERE}/run_all.sh" full > "${LAUNCH_LOG}" 2>&1 &
PID=$!
disown

{
  echo "PID=${PID}"
  echo "LOG=${LAUNCH_LOG}"
  echo "STARTED=$(date -Iseconds)"
} > "${INFO_FILE}"

echo "Pipeline launched."
echo "  PID:  ${PID}"
echo "  Log:  ${LAUNCH_LOG}"
echo "  Info: ${INFO_FILE}"
