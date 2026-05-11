#!/bin/bash
# =============================================================================
# 07_starsolo_align.sh — STARsolo alignment of 10x v3 short-read FASTQ
#
# Per GSM, we:
#   1. Build a per-mouse STARsolo CB whitelist from the authors' obs.tsv.
#      (Uses ~6,500 real cell barcodes per mouse — better than no whitelist
#      because it lets STARsolo Hamming-1 correct vs known cells.)
#   2. Run STARsolo with 10x Chromium 3' v3 chemistry:
#        R1 = 28 bp = CB(16) + UMI(12)
#        R2 = ~91 bp cDNA
#   3. Output a sorted, indexed BAM with CB / UB tags.
#
# Multi-lane handling: each GSM has 2 SRR runs (paired Element AVITI lanes).
# We pass both R2/R1 file pairs to STARsolo as comma-separated lists so it
# emits a single combined BAM per GSM.
#
# Memory: the mm10 STAR index is ~30 GB → STAR uses ~30 GB RAM during loading.
# This box has 125 GB so we run sequentially (1 GSM at a time).
#
# Resume-friendly: if {GSM}.bam already exists and is non-empty, skip.
#
# Outputs:
#   ${SR_BAM_DIR_10X}/{GSM}.bam        (sorted, CB-tagged, indexed)
#   ${SR_BAM_DIR_10X}/{GSM}.bam.bai
#   ${SR_BAM_DIR_10X}/_starsolo/{GSM}/Solo.out/  (count matrix — informational)
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/07_starsolo_align_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "07_starsolo_align.sh — STARsolo on 3 mouse 10x libraries"

activate_env
which STAR samtools >/dev/null

STAR_INDEX="${BASE_DIR}/reference/STAR_mm10_index"
if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "ERROR: STAR mm10 index missing at ${STAR_INDEX}"
    echo "Build it once via setup_references.sh"
    exit 1
fi
echo "Using STAR index: ${STAR_INDEX}"

FQ_DIR="${SR_BAM_DIR_10X}/_fastq"
SOLO_DIR="${SR_BAM_DIR_10X}/_starsolo"
WL_DIR="${SR_BAM_DIR_10X}/_whitelist"
mkdir -p "${SOLO_DIR}" "${WL_DIR}"

# -----------------------------------------------------------------------------
# Helper: build per-mouse whitelist from obs.tsv
# -----------------------------------------------------------------------------
build_whitelist_for_gsm() {
    local GSM="$1"
    local OUT="${WL_DIR}/${GSM}_whitelist.txt"
    if [[ -s "${OUT}" ]]; then
        echo "  [skip whitelist] ${OUT} exists"
        return
    fi
    local OBS
    OBS=$(ls "${BASE_DIR}/metadata/${GSM}_shortread_"*"_obs.tsv" 2>/dev/null | head -1)
    if [[ -z "${OBS}" ]]; then
        echo "ERROR: obs.tsv for ${GSM} not found in ${BASE_DIR}/metadata/"
        exit 1
    fi
    python "${THIS_DIR}/_build_starsolo_whitelist.py" \
        --obs "${OBS}" \
        --output "${OUT}"
}

# -----------------------------------------------------------------------------
# Helper: comma-join FASTQ files for a GSM (in deterministic SRR order)
# -----------------------------------------------------------------------------
fastq_list_for_gsm() {
    local GSM="$1"
    local READ_NUM="$2"  # 1 or 2
    local SRRS="${SR_10X_GSM_TO_SRRS[${GSM}]}"
    local LIST=""
    for SRR in ${SRRS}; do
        local FQ="${FQ_DIR}/${SRR}_${READ_NUM}.fastq.gz"
        if [[ ! -s "${FQ}" ]]; then
            echo "ERROR: missing FASTQ ${FQ}" >&2
            return 1
        fi
        if [[ -z "${LIST}" ]]; then
            LIST="${FQ}"
        else
            LIST="${LIST},${FQ}"
        fi
    done
    echo "${LIST}"
}

# -----------------------------------------------------------------------------
# STARsolo main loop
# -----------------------------------------------------------------------------
for GSM in "${SR_10X_GSMS[@]}"; do
    OUT_BAM="${SR_BAM_DIR_10X}/${GSM}.bam"
    if [[ -s "${OUT_BAM}" ]]; then
        echo ""
        echo "[skip] ${OUT_BAM} already exists ($(du -h "${OUT_BAM}" | cut -f1))"
        continue
    fi

    echo ""
    echo "==========================================================="
    echo "STARsolo: ${GSM}  ($(date))"
    echo "==========================================================="

    build_whitelist_for_gsm "${GSM}"
    WHITELIST="${WL_DIR}/${GSM}_whitelist.txt"

    R2_LIST=$(fastq_list_for_gsm "${GSM}" 2)
    R1_LIST=$(fastq_list_for_gsm "${GSM}" 1)
    echo "  R2 (cDNA):     ${R2_LIST}"
    echo "  R1 (CB+UMI):   ${R1_LIST}"
    echo "  Whitelist:     ${WHITELIST} ($(wc -l < "${WHITELIST}") barcodes)"

    PER_GSM_DIR="${SOLO_DIR}/${GSM}"
    mkdir -p "${PER_GSM_DIR}"

    # STARsolo flags for 10x v3 / v3.1 / v3.4 chemistry:
    #   --soloType CB_UMI_Simple             — single CB + UMI on R1
    #   --soloCBstart 1 --soloCBlen 16       — CB at R1 positions 1-16
    #   --soloUMIstart 17 --soloUMIlen 12    — UMI at R1 positions 17-28
    #   --soloFeatures Gene                  — gene-level counts
    #   --soloUMIfiltering MultiGeneUMI_CR   — filter UMIs assigned to >1 gene (CR-style)
    #   --soloUMIdedup 1MM_CR                — UMI dedup with 1-mismatch tolerance (CR-style)
    #   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts — 10x-style CB matching to whitelist
    #
    # Output:
    #   --outSAMtype BAM SortedByCoordinate
    #   --outSAMattributes NH HI AS nM CB UB GX GN
    STAR \
        --runMode alignReads \
        --runThreadN "${SR_NUM_THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${R2_LIST}" "${R1_LIST}" \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "${WHITELIST}" \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen 12 \
        --soloFeatures Gene \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM CB UB GX GN \
        --outBAMsortingThreadN 4 \
        --limitBAMsortRAM 30000000000 \
        --outTmpDir "/tmp/STAR_${GSM}_$$" \
        --outFileNamePrefix "${PER_GSM_DIR}/"

    # Move sorted BAM to the canonical location and index
    mv "${PER_GSM_DIR}/Aligned.sortedByCoord.out.bam" "${OUT_BAM}"
    samtools index -@ "${SR_NUM_THREADS}" "${OUT_BAM}"

    # Sanity check tags
    TAG_PROBE=$(samtools view "${OUT_BAM}" | head -1000 | tr '\t' '\n' | grep -oE '^[A-Z][A-Z]:' | sort -u | head -10 | tr '\n' ' ')
    echo "  Tag probe (first 1k reads): ${TAG_PROBE}"

    # Cleanup STARsolo intermediate FASTQ-derived sort files (NOT the Solo.out matrix)
    rm -f "${PER_GSM_DIR}"/_STARtmp/* 2>/dev/null || true
    rmdir "${PER_GSM_DIR}/_STARtmp" 2>/dev/null || true

    echo "  Done: ${OUT_BAM} ($(du -h "${OUT_BAM}" | cut -f1))"
done

log_step "07_starsolo_align.sh complete"
echo ""
echo "BAMs:"
ls -lh "${SR_BAM_DIR_10X}"/*.bam 2>&1 | head -10
echo ""
echo "Next: bash 08_subset_locus.sh"
