#!/bin/bash
# Run junction-PSI on every completed BAM. Concise summary table.
set -eo pipefail
THIS_DIR="$(dirname "$(readlink -f "$0")")"
source "${THIS_DIR}/00_config.sh"
activate_env

echo "Running PSI on all completed cell-line BAMs..."
echo ""
printf "%-15s %-10s %12s %5s %5s %5s %8s %12s\n" "sample" "condition" "n_at_locus" "Inc1" "Inc2" "Skip" "PSI" "Wilson_95%"
echo "----------------------------------------------------------------------------------"

for BAM in "${SR_BAM_DIR_BULK}"/*.bam; do
    [[ -s "${BAM}" ]] || continue
    SAMPLE=$(basename "${BAM}" .bam)

    # Map GSM to condition label
    case "${SAMPLE}" in
        GSM8465528) COND="NTrep1" ;;
        GSM8465529) COND="NTrep2" ;;
        GSM8465530) COND="NTrep3" ;;
        GSM8465531) COND="B15rep1" ;;
        GSM8465532) COND="B15rep2" ;;
        GSM8465533) COND="B15rep3" ;;
        GSM8465534) COND="B60rep1" ;;
        GSM8465535) COND="B60rep2" ;;
        GSM8465536) COND="B60rep3" ;;
        *) COND="?" ;;
    esac

    OUT="/tmp/${SAMPLE}_psi.tsv"
    if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
        samtools index -@ 4 "${BAM}" 2>/dev/null
    fi

    # Run calculator silently, parse output
    python "${THIS_DIR}/03_count_junctions.py" \
        --bam "${BAM}" --output "${OUT}" --sample-name "${SAMPLE}" 2>/dev/null

    # Read the row (skip header)
    read -r N_LOC INC1 INC2 SKIP PSI LO HI < <(awk -F'\t' 'NR==2 {print $2,$3,$4,$5,$6,$7,$8}' "${OUT}")
    if [[ -n "${PSI}" && "${PSI}" != "" ]]; then
        PSI_PCT=$(awk -v p="${PSI}" 'BEGIN {printf "%.3f%%", p*100}')
        LO_PCT=$(awk -v p="${LO}" 'BEGIN {printf "%.2f", p*100}')
        HI_PCT=$(awk -v p="${HI}" 'BEGIN {printf "%.2f", p*100}')
        CI="(${LO_PCT}-${HI_PCT})%"
    else
        PSI_PCT="undef"
        CI="-"
    fi
    printf "%-15s %-10s %12s %5s %5s %5s %8s %12s\n" \
        "${SAMPLE}" "${COND}" "${N_LOC}" "${INC1}" "${INC2}" "${SKIP}" "${PSI_PCT}" "${CI}"
done

echo ""
echo "Note: 95% CI = Wilson binomial CI (handles small n correctly)"
