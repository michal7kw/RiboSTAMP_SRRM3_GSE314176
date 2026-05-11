#!/bin/bash
# =============================================================================
# run_all.sh — Week 1 orchestrator
#
# Runs the full short-read junction-PSI workflow:
#   01_fetch_data.sh         — pull short-read 10x + cell-line bulk BAMs
#   04_validate_on_cell_line — validate pipeline recovers ~57% PSI
#   05_per_cell_psi_short_read — answer the progenitor question
#
# (Step 02 — locus FASTA — and step 03 — junction-PSI calculator — are
# library code called by 04 and 05, not standalone steps.)
#
# Usage:
#   bash run_all.sh                  # full run
#   bash run_all.sh validate-only    # only steps 01 + 04
#   bash run_all.sh progenitor-only  # only steps 01 + 05 (requires fetch first)
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

MODE="${1:-full}"

log_step "WEEK 1 — Short-read junction-PSI workflow (${MODE})"

case "${MODE}" in
    full)
        bash "${THIS_DIR}/01_fetch_data.sh"
        bash "${THIS_DIR}/04_validate_on_cell_line.sh"
        bash "${THIS_DIR}/05_per_cell_psi_short_read.sh"
        ;;
    validate-only)
        bash "${THIS_DIR}/01_fetch_data.sh"
        bash "${THIS_DIR}/04_validate_on_cell_line.sh"
        ;;
    progenitor-only)
        bash "${THIS_DIR}/01_fetch_data.sh"
        bash "${THIS_DIR}/05_per_cell_psi_short_read.sh"
        ;;
    *)
        echo "Usage: bash run_all.sh [full|validate-only|progenitor-only]"
        exit 1
        ;;
esac

log_step "Week 1 workflow complete"
echo ""
echo "Key outputs:"
echo "  ${SR_RESULTS_DIR}/cell_line_validation.tsv     — pipeline validation"
echo "  ${SR_RESULTS_DIR}/cell_line_validation_verdict.txt"
echo "  ${SR_RESULTS_DIR}/per_cluster_psi_short_read.tsv  — progenitor question"
echo "  ${SR_RESULTS_DIR}/progenitor_verdict.md          — written verdict"
echo ""
echo "Read progenitor_verdict.md for the answer to the PI's progenitor question."
