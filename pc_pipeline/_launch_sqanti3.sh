#!/bin/bash
# Detach setup_sqanti3.sh as a background process. See _launch_full.sh
# for why a separate launcher script is needed.
set -e

HERE="$(dirname "$(readlink -f "$0")")"
LOG=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/logs/sqanti3_setup.log
mkdir -p "$(dirname "$LOG")"

source "$HOME/miniconda3/etc/profile.d/conda.sh"

nohup bash "$HERE/setup_sqanti3.sh" > "$LOG" 2>&1 &
PID=$!
disown
echo "SQANTI3 setup launched."
echo "  PID: $PID"
echo "  Log: $LOG"
