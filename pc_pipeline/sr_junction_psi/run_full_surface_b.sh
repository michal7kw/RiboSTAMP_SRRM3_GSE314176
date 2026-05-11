#!/bin/bash
# =============================================================================
# run_full_surface_b.sh — end-to-end orchestrator for the Surface B answer
#
# Runs all 4 stages in sequence with checkpoint/resume support. Each stage
# skips work that's already done.
#
# Stage 06 — fetch 10x FASTQ from SRA (~150 GB, ~6 hours over good DSL)
# Stage 07 — STARsolo align (~30 GB RAM, ~3-4 hours per mouse, sequential)
# Stage 08 — subset to Srrm3 locus (~30 min total)
# Stage 09 — per-cell + per-cluster PSI + verdict (~5 min total)
#
# Total expected wall time: ~16-24 hours (depends on download speed).
# Total disk peak: ~250 GB during alignment, ~5 GB after cleanup.
#
# Run with `bash run_full_surface_b.sh` (or `nohup bash run_full_surface_b.sh &`
# if you want to disconnect). Resumable — if interrupted, just re-run.
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"

echo "=========================================================="
echo "Surface B — full short-read 10x progenitor analysis"
echo "Started: $(date)"
echo "=========================================================="
echo ""

bash "${THIS_DIR}/06_fetch_10x_fastq.sh"
bash "${THIS_DIR}/07_starsolo_align.sh"
bash "${THIS_DIR}/08_subset_locus.sh"
bash "${THIS_DIR}/09_per_cell_progenitor_psi.sh"

echo ""
echo "=========================================================="
echo "Surface B complete: $(date)"
echo "Verdict: $(dirname "${THIS_DIR}")/../data/sr_junction_psi/results/progenitor_verdict.md"
echo "=========================================================="
