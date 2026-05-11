#!/bin/bash
# =============================================================================
# 02_align_bulk_with_star.sh — STAR alignment for cell-line bulk samples
#
# Pulls SRA → FASTQ → STAR-aligned BAM for the 9 cell-line bulk samples.
# Output BAMs go into ${SR_BAM_DIR_BULK}/{GSM}.bam, ready for the junction-PSI
# calculator.
#
# Why this exists: SRA returns raw FASTQ for these samples (Illumina NovaSeq,
# single-end). To compute junction-PSI we need splice-aware alignment; STAR
# is the standard for splice-aware short-read RNA-seq alignment.
#
# Pre-requisites:
#   - STAR installed (`conda activate longread_iso` should have it; if not:
#     `mamba install -n longread_iso star`)
#   - mm10 STAR genome index at ${STAR_GENOME_DIR}
#   - mm10 GTF annotation at ${GENCODE_GTF}
#
# Outputs:
#   ${SR_BAM_DIR_BULK}/{GSM}.bam       — sorted, indexed STAR BAM
#   ${SR_BAM_DIR_BULK}/_star_logs/     — STAR logs
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/02_align_bulk_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "02_align_bulk_with_star.sh — fetch + align cell-line bulk samples"

# Activate conda env via helper from main 00_config.sh
activate_env

# Verify tools
which prefetch fasterq-dump STAR samtools >/dev/null

STAR_GENOME_DIR="${STAR_GENOME_DIR:-${BASE_DIR}/reference/STAR_mm10_index}"
GENCODE_GTF="${GENCODE_GTF:-${BASE_DIR}/reference/gencode.vM25.annotation.gtf}"

if [[ ! -d "${STAR_GENOME_DIR}" ]]; then
    echo "ERROR: STAR mm10 index not found at ${STAR_GENOME_DIR}"
    echo "Build it with:"
    echo "  STAR --runMode genomeGenerate \\"
    echo "       --genomeDir ${STAR_GENOME_DIR} \\"
    echo "       --genomeFastaFiles ${MM10_REF} \\"
    echo "       --sjdbGTFfile ${GENCODE_GTF} \\"
    echo "       --sjdbOverhang 99 \\"
    echo "       --runThreadN ${SR_NUM_THREADS}"
    exit 1
fi

mkdir -p "${SR_BAM_DIR_BULK}/_star_logs" "${SR_BAM_DIR_BULK}/_fastq_temp"

# -----------------------------------------------------------------------------
# Function: process_one_sample
#   $1 = GSM accession (e.g. GSM8465528)
#   $2 = SRR list (space-separated, e.g. "SRR30278197 SRR30278198")
#
# Implementation notes:
#  - FASTQs left uncompressed (we delete after alignment anyway). Saves
#    ~3-4 min per FASTQ that gzip would have spent on slow /mnt/e (9P).
#  - STAR --outTmpDir is on /tmp (ext4 native) — /mnt/e doesn't support
#    Unix FIFOs which STAR needs internally. This is the fix for the
#    "could not create FIFO file" error seen on first run.
# -----------------------------------------------------------------------------
process_one_sample() {
    local GSM="$1"
    local SRRS="$2"
    local OUT_BAM="${SR_BAM_DIR_BULK}/${GSM}.bam"

    if [[ -s "${OUT_BAM}" ]]; then
        echo "  [skip] ${OUT_BAM} already exists ($(du -h "${OUT_BAM}" | cut -f1))"
        return 0
    fi

    echo ""
    log_step "  Processing ${GSM} (${SRRS})"

    # Step 1: prefetch + extract FASTQ (NOT gzipped — saves ~3-4 min/file on /mnt/e)
    local FQ_DIR="${SR_BAM_DIR_BULK}/_fastq_temp/${GSM}"
    mkdir -p "${FQ_DIR}"
    local FQS=()
    for SRR in ${SRRS}; do
        # Either uncompressed FASTQ exists (fasterq-dump output) OR a .gz left
        # over from a previous pipeline run. Either is fine for STAR with the
        # right --readFilesCommand.
        local FQ_GZ="${FQ_DIR}/${SRR}.fastq.gz"
        local FQ_RAW="${FQ_DIR}/${SRR}.fastq"
        if [[ -s "${FQ_GZ}" ]]; then
            FQS+=("${FQ_GZ}")
            continue
        fi
        if [[ -s "${FQ_RAW}" ]]; then
            FQS+=("${FQ_RAW}")
            continue
        fi
        echo "    fetching ${SRR}..."
        prefetch --max-size 50g --output-directory "${FQ_DIR}" "${SRR}"
        local SRA="${FQ_DIR}/${SRR}/${SRR}.sra"
        fasterq-dump --threads "${SR_NUM_THREADS}" --outdir "${FQ_DIR}" "${SRA}"
        rm -rf "${FQ_DIR}/${SRR}"
        FQS+=("${FQ_RAW}")
    done

    # Step 2: STAR alignment. --outTmpDir on /tmp (ext4) to avoid 9P FIFO issue.
    # Use --readFilesCommand only if the input is gzipped.
    local FASTQ_LIST=$(IFS=,; echo "${FQS[*]}")
    local READ_CMD=""
    if [[ "${FQS[0]}" == *.gz ]]; then
        READ_CMD="--readFilesCommand zcat"
    fi
    local STAR_TMP="/tmp/STAR_${GSM}_$$"
    rm -rf "${STAR_TMP}"  # STAR errors if exists

    echo "    STAR alignment: ${FASTQ_LIST}"
    echo "    STAR tmpdir:    ${STAR_TMP}"

    # shellcheck disable=SC2086
    STAR \
        --runMode alignReads \
        --genomeDir "${STAR_GENOME_DIR}" \
        --readFilesIn "${FASTQ_LIST}" \
        ${READ_CMD} \
        --outTmpDir "${STAR_TMP}" \
        --runThreadN "${SR_NUM_THREADS}" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${SR_BAM_DIR_BULK}/_star_logs/${GSM}_" \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMattributes Standard \
        --twopassMode Basic \
        --quantMode GeneCounts \
        --limitBAMsortRAM 8000000000

    # Cleanup STAR tmp
    rm -rf "${STAR_TMP}"

    # Move the STAR BAM to the expected location
    mv "${SR_BAM_DIR_BULK}/_star_logs/${GSM}_Aligned.sortedByCoord.out.bam" "${OUT_BAM}"
    samtools index -@ "${SR_NUM_THREADS}" "${OUT_BAM}"

    # Cleanup FASTQ to save disk (~5-7 GB per uncompressed FASTQ)
    rm -rf "${FQ_DIR}"

    echo "    done: $(du -h "${OUT_BAM}" | cut -f1)"
}

# -----------------------------------------------------------------------------
# Process all 9 cell-line samples
# -----------------------------------------------------------------------------
log_step "Processing 9 cell-line bulk samples"

for GSM in "${SR_BULK_GSMS[@]}"; do
    SRRS="${SR_BULK_GSM_TO_SRRS[$GSM]}"
    if [[ -z "${SRRS}" ]]; then
        echo "WARN: no SRR list for ${GSM}, skipping"
        continue
    fi
    process_one_sample "${GSM}" "${SRRS}"
done

# Cleanup _fastq_temp if empty
rmdir "${SR_BAM_DIR_BULK}/_fastq_temp" 2>/dev/null || true

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
log_step "02_align_bulk_with_star.sh complete"
echo ""
echo "Output BAMs:"
ls -lh "${SR_BAM_DIR_BULK}"/*.bam 2>/dev/null | head -20

echo ""
echo "Next: bash 04_validate_on_cell_line.sh"
