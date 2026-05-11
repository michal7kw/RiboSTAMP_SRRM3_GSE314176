#!/bin/bash
# =============================================================================
# pc_pipeline/run_all.sh — master runner.
# Executes 01 → 04 → 04b → 05 → 07 → 08 with timestamped logs.
# Aborts on the first failed step (tee preserves stdout to log + console).
#
# Track A only:    bash run_all.sh trackA
# Full pipeline:   bash run_all.sh
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
PC_DIR="$(dirname "$(readlink -f "$0")")"

mode="${1:-full}"
ts="$(date +%Y%m%d_%H%M%S)"
runlog="${PC_LOG_DIR}/run_all_${mode}_${ts}.log"
mkdir -p "${PC_LOG_DIR}"

run_step() {
  local script="$1"
  log_step "RUN: ${script}"
  if bash "${PC_DIR}/${script}" 2>&1 | tee -a "${runlog}"; then
    echo "OK: ${script}" | tee -a "${runlog}"
  else
    echo "FAIL: ${script}. See ${runlog}." | tee -a "${runlog}" >&2
    exit 1
  fi
}

case "${mode}" in
  trackA)
    run_step 01_fetch_processed.sh
    log_step "Track A complete. Now run local R scripts (local/01, 05, 06)."
    ;;
  full)
    # References (mm10 + Gencode + liftOver chains) are a one-time prereq
    # for steps 04 onward. setup_references.sh is idempotent — skips files
    # that already exist — so it's safe to call at the head of every run.
    run_step setup_references.sh
    run_step 01_fetch_processed.sh
    run_step 02_fetch_pacbio.sh
    run_step 03_skera_split.sh
    run_step 04_align_pbmm2.sh
    run_step 04b_verify_anchor_coverage.sh
    run_step 05_isoquant_targeted.sh
    run_step 07_sqanti3_qc.sh
    run_step 08_liftover_mm10_to_mm39.sh
    log_step "Full pipeline complete. Now run local R scripts (local/01-07)."
    log_step "NB: SQANTI3 (step 07) requires its own env. Run pc_pipeline/setup_sqanti3.sh once."
    ;;
  *)
    echo "Usage: $0 [trackA|full]" >&2
    exit 1
    ;;
esac

echo ""
echo "Run log: ${runlog}"
