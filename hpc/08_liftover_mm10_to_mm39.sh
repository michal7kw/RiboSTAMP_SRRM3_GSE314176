#!/bin/bash
#SBATCH --job-name=ribo_08_liftover
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/08_liftover_%j.err
#SBATCH --output=./logs/08_liftover_%j.out

# =============================================================================
# liftOver the novel-isoform GTF from mm10 to mm39, then cross-reference
# against the lab's short-read novel-exon anchor at:
#     chr5:135,898,574-135,898,652 (mm39, 79 bp, negative strand)
#
# anchor source: ../90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
SQANTI_DIR="${BASE_DIR}/data/sqanti3"
OUT_DIR="${BASE_DIR}/data/liftover"
mkdir -p "${OUT_DIR}" "$(pwd)/logs"

CHAIN="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10ToMm39.over.chain.gz"
# If chain file is not already present, provision it once:
if [ ! -f "${CHAIN}" ]; then
  echo "Chain file missing — attempting auto-download from UCSC goldenPath."
  mkdir -p "$(dirname "${CHAIN}")"
  wget -O "${CHAIN}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
fi

NOVEL_GTF="${SQANTI_DIR}/ribostamp_srrm3_corrected.gtf"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

echo "=== Step 1: GTF -> genePred -> BED12 (liftOver wants BED) ==="
gtfToGenePred "${NOVEL_GTF}" "${OUT_DIR}/novel_mm10.gp"
genePredToBed "${OUT_DIR}/novel_mm10.gp" "${OUT_DIR}/novel_mm10.bed"

echo "=== Step 2: liftOver mm10 -> mm39 ==="
liftOver \
    "${OUT_DIR}/novel_mm10.bed" \
    "${CHAIN}" \
    "${OUT_DIR}/novel_mm39.bed" \
    "${OUT_DIR}/novel_mm39.unmapped.bed"

echo "=== Step 3: cross-ref against short-read anchor ==="
# Anchor: chr5:135898574-135898652 (mm39), 79 bp, negative strand.
printf "chr5\t135898574\t135898652\tshortread_anchor_79bp\t.\t-\n" \
    > "${OUT_DIR}/anchor_mm39.bed"

echo "--- direct overlap ---"
bedtools intersect \
    -a "${OUT_DIR}/novel_mm39.bed" \
    -b "${OUT_DIR}/anchor_mm39.bed" \
    -wa -wb \
    > "${OUT_DIR}/overlap_with_anchor.bed" || true

if [ -s "${OUT_DIR}/overlap_with_anchor.bed" ]; then
  echo "MATCH: novel isoform(s) overlap the short-read 79-bp anchor."
  cat "${OUT_DIR}/overlap_with_anchor.bed"
else
  echo "NO DIRECT OVERLAP. Computing nearest-feature distance..."
  bedtools sort -i "${OUT_DIR}/novel_mm39.bed"    > "${OUT_DIR}/novel_mm39.sorted.bed"
  bedtools sort -i "${OUT_DIR}/anchor_mm39.bed"   > "${OUT_DIR}/anchor_mm39.sorted.bed"
  bedtools closest \
      -a "${OUT_DIR}/anchor_mm39.sorted.bed" \
      -b "${OUT_DIR}/novel_mm39.sorted.bed" \
      -d \
      > "${OUT_DIR}/nearest_novel_to_anchor.bed"
  echo "Nearest novel isoforms (last column = distance in bp):"
  cat "${OUT_DIR}/nearest_novel_to_anchor.bed"
fi

echo "=== Done: $(date) ==="
