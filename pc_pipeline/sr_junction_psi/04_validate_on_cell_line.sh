#!/bin/bash
# =============================================================================
# 04_validate_on_cell_line.sh — pipeline validation on cell-line bulk data
#
# Runs the junction-PSI calculator on the 9 cell-line bulk samples (NT/B15/B60).
# Expected outcome: NT samples (rep 1-3) should give PSI ≈ 50-60% (matching
# the lab's bulk anchor study), confirming the pipeline can detect inclusion
# when present. If we get ~0%, the pipeline has a silent bug.
#
# Outputs:
#   ${SR_RESULTS_DIR}/cell_line_validation.tsv  — per-sample PSI table
#   ${SR_RESULTS_DIR}/cell_line_validation_verdict.txt — PASS / FAIL banner
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/04_validate_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "04_validate_on_cell_line.sh — checking pipeline against ground truth"

activate_env

# Pre-flight: do we have the bulk BAMs?
N_BAMS=$(ls "${SR_BAM_DIR_BULK}"/*.bam 2>/dev/null | wc -l)
if [[ "${N_BAMS}" -eq 0 ]]; then
    echo "ERROR: no cell-line BAMs in ${SR_BAM_DIR_BULK}"
    echo "Run 01_fetch_data.sh first (and confirm SR_BULK_SRRS is populated)"
    exit 1
fi
echo "Found ${N_BAMS} cell-line BAMs"

# Run junction-PSI on each
PER_SAMPLE_TSV="${SR_RESULTS_DIR}/cell_line_validation.tsv"
> "${PER_SAMPLE_TSV}"
HEADER_DONE=0

for BAM in "${SR_BAM_DIR_BULK}"/*.bam; do
    SAMPLE=$(basename "${BAM}" .bam)
    OUT="${SR_PER_CELL_DIR}/${SAMPLE}_bulk_psi.tsv"
    echo ""
    echo "Processing ${SAMPLE}..."

    python "${THIS_DIR}/03_count_junctions.py" \
        --bam "${BAM}" \
        --output "${OUT}" \
        --sample-name "${SAMPLE}"

    # Append to combined TSV (keep header from first sample only)
    if [[ "${HEADER_DONE}" -eq 0 ]]; then
        cat "${OUT}" >> "${PER_SAMPLE_TSV}"
        HEADER_DONE=1
    else
        tail -n +2 "${OUT}" >> "${PER_SAMPLE_TSV}"
    fi
done

# -----------------------------------------------------------------------------
# Verdict — did NT samples recover ~57% PSI?
# -----------------------------------------------------------------------------
log_step "Verdict — comparing NT samples to lab's bulk anchor of 57%"

VERDICT="${SR_RESULTS_DIR}/cell_line_validation_verdict.txt"
python <<EOF | tee "${VERDICT}"
import pandas as pd
df = pd.read_csv("${PER_SAMPLE_TSV}", sep="\t")
print()
print("Per-sample bulk PSI:")
print(df[["sample", "n_inc1", "n_inc2", "n_skip", "psi", "psi_lo95", "psi_hi95"]].to_string(index=False))

# Identify NT samples (the 'NT' samples are the WT condition expected to have ~57% PSI)
nt_mask = df["sample"].str.contains("NT", case=False, na=False)
if nt_mask.sum() == 0:
    print()
    print("WARN: no NT samples found in BAM list (sample names should contain 'NT' for the WT condition)")
    print("VERDICT: INCONCLUSIVE — could not identify validation samples")
else:
    nt_psi = df.loc[nt_mask, "psi"].dropna()
    if len(nt_psi) == 0:
        print()
        print("VERDICT: FAIL — all NT samples have undefined PSI (no informative reads at the locus)")
    else:
        mean_nt = nt_psi.mean() * 100
        print()
        print(f"Mean PSI in NT samples: {mean_nt:.1f}% (n={len(nt_psi)} samples)")
        print(f"Lab's bulk anchor target: ~57% (Parental WT)")
        if 30 <= mean_nt <= 75:
            print(f"VERDICT: PASS — pipeline recovers cassette inclusion in cell-line samples")
            print(f"  → 0% PSI in adult hippocampus is real biology, not a pipeline artifact")
        elif mean_nt < 30:
            print(f"VERDICT: FAIL — NT mean PSI ({mean_nt:.1f}%) is far below expected 57%")
            print(f"  → pipeline may have a silent bug; investigate before trusting hippocampus 0%")
        else:
            print(f"VERDICT: WARN — NT mean PSI ({mean_nt:.1f}%) is above 75% (unexpected high)")
            print(f"  → pipeline may be over-counting inclusion; investigate")
EOF

log_step "04_validate_on_cell_line.sh complete"
echo ""
echo "Outputs:"
echo "  Per-sample table: ${PER_SAMPLE_TSV}"
echo "  Verdict:          ${VERDICT}"
