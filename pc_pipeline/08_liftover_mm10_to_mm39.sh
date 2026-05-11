#!/bin/bash
# =============================================================================
# pc_pipeline/08 — liftOver novel-isoform GTF mm10 → mm39, cross-ref anchor.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

NOVEL_GTF="${SQANTI_DIR}/ribostamp_srrm3_corrected.gtf"
if [ ! -f "${NOVEL_GTF}" ]; then
  echo "ERROR: ${NOVEL_GTF} not found; run 07_sqanti3_qc.sh first." >&2
  exit 1
fi

if [ ! -f "${CHAIN_MM10_TO_MM39}" ]; then
  log_step "Fetching mm10ToMm39 chain (one-time)"
  mkdir -p "${CHAIN_DIR}"
  wget -O "${CHAIN_MM10_TO_MM39}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
fi

log_step "Step 1 — GTF -> genePred -> BED12"
gtfToGenePred "${NOVEL_GTF}" "${LIFT_DIR}/novel_mm10.gp"
genePredToBed "${LIFT_DIR}/novel_mm10.gp" "${LIFT_DIR}/novel_mm10.bed"

log_step "Step 2 — liftOver mm10 -> mm39"
liftOver \
    "${LIFT_DIR}/novel_mm10.bed" \
    "${CHAIN_MM10_TO_MM39}" \
    "${LIFT_DIR}/novel_mm39.bed" \
    "${LIFT_DIR}/novel_mm39.unmapped.bed"

log_step "Step 3 — cross-ref against short-read anchor"
printf "chr5\t135898574\t135898652\tshortread_anchor_79bp\t.\t-\n" > "${LIFT_DIR}/anchor_mm39.bed"

bedtools intersect \
    -a "${LIFT_DIR}/novel_mm39.bed" \
    -b "${LIFT_DIR}/anchor_mm39.bed" \
    -wa -wb > "${LIFT_DIR}/overlap_with_anchor.bed" || true

if [ -s "${LIFT_DIR}/overlap_with_anchor.bed" ]; then
  echo "MATCH: novel isoform(s) overlap the short-read 79-bp anchor."
  cat "${LIFT_DIR}/overlap_with_anchor.bed"
else
  echo "NO DIRECT OVERLAP. Computing nearest-feature distance..."
  bedtools sort -i "${LIFT_DIR}/novel_mm39.bed"  > "${LIFT_DIR}/novel_mm39.sorted.bed"
  bedtools sort -i "${LIFT_DIR}/anchor_mm39.bed" > "${LIFT_DIR}/anchor_mm39.sorted.bed"
  bedtools closest \
      -a "${LIFT_DIR}/anchor_mm39.sorted.bed" \
      -b "${LIFT_DIR}/novel_mm39.sorted.bed" \
      -d > "${LIFT_DIR}/nearest_novel_to_anchor.bed"
  echo "Nearest novel isoforms (last column = distance bp):"
  cat "${LIFT_DIR}/nearest_novel_to_anchor.bed"
fi

log_step "Done"
