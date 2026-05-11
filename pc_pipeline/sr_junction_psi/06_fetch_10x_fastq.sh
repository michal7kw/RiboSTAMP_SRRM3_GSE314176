#!/bin/bash
# =============================================================================
# 06_fetch_10x_fastq.sh — download 10x short-read FASTQ from SRA
#
# The GEO/SRA submission for GSE314176 deposits the 10x scRNA-seq runs as
# RAW FASTQ (no aligned BAM was uploaded). We must therefore download the
# paired R1/R2 FASTQ for each SRR and align it ourselves with STARsolo.
#
# Per GSM (mouse): 2 SRR runs (paired Element AVITI lanes).
# Per SRR: ~39 GB SRA archive → ~80 GB paired FASTQ (gzipped ~25 GB).
# Total disk for FASTQ: ~3 × 2 × 25 GB = ~150 GB across all 3 mice.
#
# Resume-friendly: if the FASTQ pair for a SRR already exists on disk and is
# non-empty, the script skips that SRR.
#
# Output layout:
#   ${SR_BAM_DIR_10X}/_fastq/{SRR}_1.fastq.gz   (R1: CB + UMI)
#   ${SR_BAM_DIR_10X}/_fastq/{SRR}_2.fastq.gz   (R2: cDNA)
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/06_fetch_10x_fastq_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "06_fetch_10x_fastq.sh — pulling 10x short-read FASTQ from SRA"

activate_env
which prefetch fasterq-dump >/dev/null

FQ_DIR="${SR_BAM_DIR_10X}/_fastq"
SRA_DIR="${SR_BAM_DIR_10X}/_prefetch"
mkdir -p "${FQ_DIR}" "${SRA_DIR}"

# Pre-flight: free disk
FREE_GB=$(df -BG "${BASE_DIR}" | awk 'NR==2 {print $4}' | tr -d 'G')
echo "Free disk on ${BASE_DIR}: ${FREE_GB} GB"
if [[ "${FREE_GB}" -lt 200 ]]; then
    echo "WARN: <200 GB free; FASTQ download needs ~150 GB. Continuing anyway."
fi

# -----------------------------------------------------------------------------
# Function: fetch_one_srr
#   $1 = SRR accession
# -----------------------------------------------------------------------------
fetch_one_srr() {
    local SRR="$1"
    local R1="${FQ_DIR}/${SRR}_1.fastq.gz"
    local R2="${FQ_DIR}/${SRR}_2.fastq.gz"

    if [[ -s "${R1}" && -s "${R2}" ]]; then
        echo "  [skip] ${SRR}: FASTQ pair already exists ($(du -h "${R1}" | cut -f1) + $(du -h "${R2}" | cut -f1))"
        return 0
    fi

    echo ""
    echo "  ----- ${SRR} -----"
    echo "  $(date '+%H:%M:%S') prefetch ${SRR}"

    local SRA_FILE="${SRA_DIR}/${SRR}/${SRR}.sra"
    if [[ ! -s "${SRA_FILE}" ]]; then
        prefetch --max-size 100g --output-directory "${SRA_DIR}" "${SRR}"
    else
        echo "  [skip prefetch] SRA archive already exists"
    fi

    echo "  $(date '+%H:%M:%S') fasterq-dump ${SRR} (paired, threaded)"
    # --split-files: write _1.fastq + _2.fastq for paired data
    # -e: threads (use most of them — fasterq-dump is I/O bound but parallel)
    # --skip-technical: drop the 10x sample-index read if present
    # -O: output dir
    fasterq-dump \
        --split-files \
        --skip-technical \
        --threads "${SR_NUM_THREADS}" \
        --temp "${SRA_DIR}/_tmp_${SRR}" \
        --outdir "${FQ_DIR}" \
        "${SRA_FILE}"

    echo "  $(date '+%H:%M:%S') gzipping FASTQ"
    # gzip in parallel (foreground so we don't oversubscribe)
    pigz -p "${SR_NUM_THREADS}" "${FQ_DIR}/${SRR}_1.fastq" "${FQ_DIR}/${SRR}_2.fastq"

    # Cleanup intermediate
    rm -rf "${SRA_DIR}/${SRR}" "${SRA_DIR}/_tmp_${SRR}"
    echo "  $(date '+%H:%M:%S') done ${SRR}: $(du -h "${R1}" | cut -f1) + $(du -h "${R2}" | cut -f1)"
}

# -----------------------------------------------------------------------------
# Loop: fetch all SRRs from the GSM-to-SRRs map
# -----------------------------------------------------------------------------
for GSM in "${SR_10X_GSMS[@]}"; do
    SRRS="${SR_10X_GSM_TO_SRRS[${GSM}]}"
    echo ""
    echo "=========================================================="
    echo "GSM ${GSM}: SRRs = ${SRRS}"
    echo "=========================================================="
    for SRR in ${SRRS}; do
        fetch_one_srr "${SRR}"
    done
done

log_step "06_fetch_10x_fastq.sh complete"
echo ""
ls -lh "${FQ_DIR}"/ 2>&1 | head -20
echo ""
echo "Total FASTQ size: $(du -sh "${FQ_DIR}" | cut -f1)"
echo ""
echo "Next: bash 07_starsolo_align.sh"
