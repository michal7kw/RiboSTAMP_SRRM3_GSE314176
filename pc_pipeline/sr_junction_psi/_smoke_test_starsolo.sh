#!/bin/bash
# =============================================================================
# _smoke_test_starsolo.sh — verify STARsolo config before full-scale runs
#
# Takes the first 1M paired reads from the smallest already-downloaded FASTQ,
# runs STARsolo against the per-mouse whitelist, and reports:
#   - mapping rate
#   - presence of CB/UB tags in BAM
#   - number of reads at the Srrm3 locus
#   - the per-cell BAM CB tag value distribution
#
# Catches config errors (wrong read structure, wrong whitelist orientation,
# missing flags) BEFORE we waste hours on the full-scale run.
#
# Run after at least one SRR's FASTQ pair has finished pigz'ing.
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

activate_env

SMOKE_DIR="${SR_BAM_DIR_10X}/_smoke_test"
mkdir -p "${SMOKE_DIR}"

FQ_DIR="${SR_BAM_DIR_10X}/_fastq"
FIRST_R1=$(ls -S "${FQ_DIR}"/*_1.fastq.gz 2>/dev/null | tail -1)
if [[ -z "${FIRST_R1}" ]]; then
    echo "ERROR: no FASTQ pairs ready in ${FQ_DIR}"
    exit 1
fi
SRR=$(basename "${FIRST_R1}" _1.fastq.gz)
R2="${FQ_DIR}/${SRR}_2.fastq.gz"
echo "Smoke-testing on SRR=${SRR}"
echo "  R1: ${FIRST_R1}"
echo "  R2: ${R2}"

# Find which GSM this SRR belongs to (from config map)
GSM=""
for G in "${SR_10X_GSMS[@]}"; do
    SRRS="${SR_10X_GSM_TO_SRRS[${G}]}"
    if [[ " ${SRRS} " == *" ${SRR} "* ]]; then
        GSM="${G}"
        break
    fi
done
if [[ -z "${GSM}" ]]; then
    echo "ERROR: could not match SRR ${SRR} to a GSM"
    exit 1
fi
echo "  GSM = ${GSM}"

# Build whitelist for this mouse
WL="${SMOKE_DIR}/${GSM}_whitelist.txt"
if [[ ! -s "${WL}" ]]; then
    OBS=$(ls "${BASE_DIR}/metadata/${GSM}_shortread_"*"_obs.tsv" 2>/dev/null | head -1)
    python "${THIS_DIR}/_build_starsolo_whitelist.py" --obs "${OBS}" --output "${WL}"
fi
echo "  Whitelist: $(wc -l < "${WL}") barcodes"

# Subsample first 1M reads (4M lines per FASTQ).
# `head` closes its stdin once it has the first N lines, sending SIGPIPE to
# zcat upstream. With `set -o pipefail` that propagates as a non-zero exit
# even though we got the data we wanted. Disable pipefail just for these.
R1_SMALL="${SMOKE_DIR}/${SRR}_1.small.fastq.gz"
R2_SMALL="${SMOKE_DIR}/${SRR}_2.small.fastq.gz"
if [[ ! -s "${R1_SMALL}" ]]; then
    echo "Subsampling first 1M reads..."
    set +o pipefail
    zcat "${FIRST_R1}" | head -n 4000000 | pigz -p 8 > "${R1_SMALL}"
    zcat "${R2}"        | head -n 4000000 | pigz -p 8 > "${R2_SMALL}"
    set -o pipefail
fi
echo "  R1 subset: $(du -h "${R1_SMALL}" | cut -f1)"
echo "  R2 subset: $(du -h "${R2_SMALL}" | cut -f1)"

# Inspect read length structure (disable pipefail — head closes upstream zcat → SIGPIPE).
echo ""
echo "--- Read length sanity ---"
set +o pipefail
echo "R1 (expect 28bp = CB16 + UMI12):"
zcat "${R1_SMALL}" | head -8 | awk 'NR%4==2 {print "  len=" length($0) "  seq=" substr($0,1,30) "..."}'
echo "R2 (expect ~91bp cDNA):"
zcat "${R2_SMALL}" | head -8 | awk 'NR%4==2 {print "  len=" length($0) "  seq=" substr($0,1,30) "..."}'
set -o pipefail

# Run STARsolo on the subset
OUT_DIR="${SMOKE_DIR}/${GSM}_${SRR}_smoke"
mkdir -p "${OUT_DIR}"
STAR_INDEX="${BASE_DIR}/reference/STAR_mm10_index"

echo ""
echo "--- Running STARsolo on 1M-read subset ---"
STAR \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${R2_SMALL}" "${R1_SMALL}" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WL}" \
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
    --outTmpDir "/tmp/STAR_smoke_$$" \
    --outFileNamePrefix "${OUT_DIR}/"

BAM="${OUT_DIR}/Aligned.sortedByCoord.out.bam"
samtools index -@ 8 "${BAM}"

# Verify BAM
echo ""
echo "--- BAM verification ---"
echo "BAM size: $(du -h "${BAM}" | cut -f1)"
echo "Total reads in BAM: $(samtools view -c "${BAM}")"
echo "Reads with CB tag: $(samtools view "${BAM}" | grep -c 'CB:Z:' || true)"
echo "Reads at Srrm3 locus chr5:135680000-135960000:"
samtools view -c "${BAM}" chr5:135680000-135960000

echo ""
echo "--- Solo.out summary ---"
ls -la "${OUT_DIR}/Solo.out/" 2>&1 | head -10
cat "${OUT_DIR}/Solo.out/Gene/Summary.csv" 2>/dev/null || echo "(no Summary.csv yet)"

echo ""
echo "--- STAR Log.final.out (mapping stats) ---"
cat "${OUT_DIR}/Log.final.out" | head -40

echo ""
echo "--- Tag distribution in first 1000 reads at locus ---"
samtools view "${BAM}" chr5:135680000-135960000 2>/dev/null | head -1000 | tr '\t' '\n' | grep -oE '^[A-Z][A-Z]:' | sort | uniq -c | sort -rn

echo ""
echo "Smoke test PASSED if:"
echo "  1. R1 length is 28, R2 length is ~91"
echo "  2. Mapping rate >= 60% in Log.final.out"
echo "  3. CB:Z: tags present in BAM"
echo "  4. Some reads at the Srrm3 locus"
