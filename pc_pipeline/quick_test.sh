#!/bin/bash
# =============================================================================
# pc_pipeline/quick_test.sh — validate the pipeline end-to-end on the smallest
# single SRR (SRR36480452, 115 GB). Total wall-clock ~6-10 hours on this PC.
#
# Use this BEFORE committing to the full 374 GB / multi-day run. It catches
# tooling problems (wrong skera adapter, missing reference, conda env issues)
# at the cheapest possible scale.
#
# What it does NOT do: produce statistically meaningful Fisher results — only
# 1 of 3 biological replicates means the progenitor enrichment test is
# underpowered. Treat the output as "the pipeline works; novel cassette is /
# is not detectable in this one replicate".
# =============================================================================

set -euo pipefail

# Override SRR_LIST to just the smallest run
export SRR_LIST="SRR36480452"

PC_DIR="$(dirname "$(readlink -f "$0")")"
log_step() { printf "\n=== [%s] %s ===\n" "$(date '+%H:%M:%S')" "$*"; }

log_step "Quick-test mode: only ${SRR_LIST}"
bash "${PC_DIR}/run_all.sh" full
