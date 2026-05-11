#!/bin/bash
# Quick PSI sanity check on the first cell-line BAM as it becomes available.
# This is a preview of what 04_validate_on_cell_line.sh will do at the end.
set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"
activate_env

echo "=== Quick junction-PSI sanity check ==="
echo "Looking for first available cell-line BAM..."

for BAM in "${SR_BAM_DIR_BULK}"/*.bam; do
    [[ -s "${BAM}" ]] || continue
    SAMPLE=$(basename "${BAM}" .bam)
    OUT="/tmp/${SAMPLE}_quickpsi.tsv"
    echo ""
    echo ">>> Sample: ${SAMPLE}"
    echo ">>> BAM:    ${BAM} ($(du -h "${BAM}" | cut -f1))"

    if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
        echo ">>> No index found; building..."
        samtools index -@ 8 "${BAM}"
    fi

    python "${THIS_DIR}/03_count_junctions.py" \
        --bam "${BAM}" \
        --output "${OUT}" \
        --sample-name "${SAMPLE}"

    echo ""
    cat "${OUT}"
    break  # Only first one for the quick check
done

echo ""
echo "=== Interpretation ==="
echo "Lab anchor study expects PSI ≈ 57% in NT samples (Parental/WT)."
echo "GSM8465528 = NTrep1 (the WT condition), so we expect ~57% PSI."
echo "If we recover 30-75%, the pipeline is validated; <30% suggests a bug."
