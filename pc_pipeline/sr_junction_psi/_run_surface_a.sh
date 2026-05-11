#!/bin/bash
# =============================================================================
# _run_surface_a.sh — Surface A end-to-end (cell-line validation)
#
# Builds STAR index (one-time, ~30-60 min), downloads 9 cell-line bulk SRA
# archives, aligns with STAR, runs junction-PSI calculator, emits verdict.
#
# Total expected wall time: 2-4 hours (depends on bandwidth, CPU).
# =============================================================================
set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/surface_a_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "Surface A — cell-line validation, full pipeline"

activate_env

STAR_GENOME_DIR="${STAR_GENOME_DIR:-${BASE_DIR}/reference/STAR_mm10_index}"
GENCODE_GTF="${GENCODE_GTF:-${BASE_DIR}/reference/mm10/gencode.vM25.annotation.gtf}"
MM10_FA="${BASE_DIR}/reference/mm10/mm10.fa"

# -----------------------------------------------------------------------------
# Step 1: Build STAR index if missing (one-time)
# -----------------------------------------------------------------------------
if [[ ! -f "${STAR_GENOME_DIR}/Genome" ]]; then
    log_step "Building STAR mm10 index (one-time, ~30-60 min)"
    mkdir -p "${STAR_GENOME_DIR}"
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_GENOME_DIR}" \
        --genomeFastaFiles "${MM10_FA}" \
        --sjdbGTFfile "${GENCODE_GTF}" \
        --sjdbOverhang 99 \
        --runThreadN "${SR_NUM_THREADS}" \
        --genomeSAsparseD 2          # smaller index, slightly slower
    echo "STAR index built at ${STAR_GENOME_DIR}"
else
    log_step "STAR index found at ${STAR_GENOME_DIR}, reusing"
fi

# -----------------------------------------------------------------------------
# Step 2: download + align 9 cell-line bulk samples
# -----------------------------------------------------------------------------
log_step "Step 2: download + align cell-line bulk samples"
bash "${THIS_DIR}/02_align_bulk_with_star.sh"

# -----------------------------------------------------------------------------
# Step 3: run junction-PSI calculator on each, emit verdict
# -----------------------------------------------------------------------------
log_step "Step 3: junction-PSI counting + verdict"
bash "${THIS_DIR}/04_validate_on_cell_line.sh"

log_step "Surface A complete"
echo ""
echo "Read the verdict at: ${SR_RESULTS_DIR}/cell_line_validation_verdict.txt"
