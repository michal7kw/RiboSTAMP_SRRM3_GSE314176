#!/bin/bash
# =============================================================================
# 01_fetch_data.sh
# Pull short-read 10x BAMs (3 mice, ~50-80 GB each) + cell-line bulk RNA-seq
# (9 small samples, ~3-5 GB each) from SRA.
#
# Two surfaces:
#   - "10x" — short-read scRNA-seq, contains OPC + iDG-Immature for the
#     progenitor question
#   - "bulk" — cell-line samples, used to validate the pipeline recovers
#     ~57% PSI at the cassette
#
# Outputs:
#   - {SR_BAM_DIR_10X}/{GSM}.bam    (10x BAMs with CB/UB tags)
#   - {SR_BAM_DIR_BULK}/{GSM}.bam   (cell-line bulk BAMs)
#
# IMPORTANT: 10x BAMs are ~50-80 GB EACH. Total disk: 150-250 GB for 3 mice
# plus 30-45 GB for the 9 cell-line samples. Verify free space before running.
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/01_fetch_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "01_fetch_data.sh — pulling short-read BAMs from SRA"

# Pre-flight: free disk
FREE_GB=$(df -BG "${BASE_DIR}" | awk 'NR==2 {print $4}' | tr -d 'G')
echo "Free disk on ${BASE_DIR}: ${FREE_GB} GB"
if [[ "${FREE_GB}" -lt 200 ]]; then
    echo "WARN: <200 GB free; 10x BAMs are large. Consider freeing space first."
fi

activate_env
which prefetch fasterq-dump sam-dump samtools >/dev/null

# -----------------------------------------------------------------------------
# Function: fetch_one_sra
#   $1 = SRR accession
#   $2 = output BAM path
#   $3 = data_kind ("10x" or "bulk")
# -----------------------------------------------------------------------------
fetch_one_sra() {
    local SRR="$1"
    local OUT_BAM="$2"
    local KIND="$3"

    if [[ -s "${OUT_BAM}" ]]; then
        echo "  [skip] ${OUT_BAM} already exists ($(du -h "${OUT_BAM}" | cut -f1))"
        return 0
    fi

    echo "  fetching ${SRR} → ${OUT_BAM} (${KIND})"

    # Step 1: prefetch the SRA archive
    local SRA_DIR
    SRA_DIR="$(dirname "${OUT_BAM}")/_prefetch"
    mkdir -p "${SRA_DIR}"
    prefetch --max-size 100g --output-directory "${SRA_DIR}" "${SRR}"

    # Step 2: convert to BAM
    # For 10x scRNA-seq: SRA has the CB/UB tags preserved in the BAM. Use sam-dump.
    # For cell-line bulk: same — sam-dump preserves alignment data.
    local SRA_FILE="${SRA_DIR}/${SRR}/${SRR}.sra"
    if [[ ! -f "${SRA_FILE}" ]]; then
        echo "ERROR: prefetch didn't produce ${SRA_FILE}"
        return 1
    fi

    sam-dump --aligned --output-file - "${SRA_FILE}" \
        | samtools view -bS -@ "${SR_NUM_THREADS}" - \
        | samtools sort -@ "${SR_NUM_THREADS}" -o "${OUT_BAM}" -
    samtools index -@ "${SR_NUM_THREADS}" "${OUT_BAM}"

    # Cleanup the SRA archive (we don't need it anymore)
    rm -rf "${SRA_DIR}/${SRR}"

    echo "  done: $(du -h "${OUT_BAM}" | cut -f1)"
}

# -----------------------------------------------------------------------------
# Fetch the 10x short-read scRNA-seq BAMs (3 mice)
# -----------------------------------------------------------------------------
log_step "10x short-read scRNA-seq (3 mice — has OPC + iDG-Immature)"
for i in "${!SR_10X_SRRS[@]}"; do
    SRR="${SR_10X_SRRS[$i]}"
    GSM="${SR_10X_GSMS[$i]}"
    fetch_one_sra "${SRR}" "${SR_BAM_DIR_10X}/${GSM}.bam" "10x"
done

# -----------------------------------------------------------------------------
# Fetch the cell-line bulk RNA-seq BAMs (9 samples — validation)
# -----------------------------------------------------------------------------
log_step "Cell-line bulk Ribo-STAMP samples (9 — validation set, expect ~57% PSI in NT)"
echo "  TODO: SR_BULK_SRRS array is not yet populated in 00_config.sh."
echo "  To populate: visit https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE275091"
echo "  and pull the SRR-to-GSM mapping for the 9 samples (48dox-RPS2-NT/B15/B60)."
echo "  Or: query SRA programmatically with edirect:"
echo "    esearch -db sra -query 'PRJNA1390205 AND 48dox' | efetch -format runinfo"
echo ""
echo "  Stub fetch loop (will run once SR_BULK_SRRS is populated):"
echo ""
# Check if SR_BULK_SRRS array exists and is non-empty (uncomment when populated)
# for i in "${!SR_BULK_SRRS[@]}"; do
#     SRR="${SR_BULK_SRRS[$i]}"
#     GSM="${SR_BULK_GSMS[$i]}"
#     fetch_one_sra "${SRR}" "${SR_BAM_DIR_BULK}/${GSM}.bam" "bulk"
# done

# -----------------------------------------------------------------------------
# Sanity check: confirm BAMs have the expected tags
# -----------------------------------------------------------------------------
log_step "Sanity check — BAM tag presence"

if [[ -d "${SR_BAM_DIR_10X}" && -n "$(ls -A "${SR_BAM_DIR_10X}"/*.bam 2>/dev/null)" ]]; then
    echo "10x BAMs — expected tags: CB (cell barcode), UB (UMI), GX/GN (gene)"
    set +e
    for BAM in "${SR_BAM_DIR_10X}"/*.bam; do
        TAGS=$(samtools view "${BAM}" | head -100 | tr '\t' '\n' | grep -oE '^[A-Z][A-Z]:' | sort -u | head -10 | tr '\n' ' ')
        echo "  ${BAM}: ${TAGS}"
    done
    set -e
fi

log_step "01_fetch_data.sh complete"
echo ""
echo "Outputs:"
echo "  10x BAMs:  ${SR_BAM_DIR_10X}/*.bam  ($(ls "${SR_BAM_DIR_10X}"/*.bam 2>/dev/null | wc -l) files)"
echo "  Bulk BAMs: ${SR_BAM_DIR_BULK}/*.bam  ($(ls "${SR_BAM_DIR_BULK}"/*.bam 2>/dev/null | wc -l) files)"
echo ""
echo "Next: bash 03_count_junctions.py (after populating SR_BULK_SRRS)"
