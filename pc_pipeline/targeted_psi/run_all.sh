#!/bin/bash
# =============================================================================
# Targeted PSI master runner.
# Runs in parallel with the main pipeline — uses TGT_NUM_THREADS (default 6)
# instead of NUM_THREADS (16) so we coexist gracefully on the i5-14600K.
# Outputs go to ${TGT_DIR}, separate from anything the main pipeline writes.
# =============================================================================

set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "${HERE}/00_config.sh"

ts="$(date +%Y%m%d_%H%M%S)"
runlog="${TGT_LOG_DIR}/run_all_${ts}.log"

run_step() {
  local script="$1"
  log_step "RUN: ${script}"
  if bash "${HERE}/${script}" 2>&1 | tee -a "${runlog}"; then
    echo "OK: ${script}" | tee -a "${runlog}"
  else
    echo "FAIL: ${script}. See ${runlog}." | tee -a "${runlog}" >&2
    exit 1
  fi
}

run_python_step() {
  local script="$1"
  log_step "RUN (python): ${script}"
  if python "${HERE}/${script}" 2>&1 | tee -a "${runlog}"; then
    echo "OK: ${script}" | tee -a "${runlog}"
  else
    echo "FAIL: ${script}. See ${runlog}." | tee -a "${runlog}" >&2
    exit 1
  fi
}

# Make sure conda env (longread_iso) is active for python scripts too
activate_env

run_step       01_extract_locus.sh
run_step       02_align_to_locus.sh
run_python_step 03_classify_and_barcode.py
run_python_step 04_per_cluster_psi.py

log_step "Targeted PSI complete."
echo "Run log: ${runlog}"
echo "Results:"
ls -lh "${TGT_RESULTS_DIR}/"
