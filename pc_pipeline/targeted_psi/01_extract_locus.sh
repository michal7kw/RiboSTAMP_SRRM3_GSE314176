#!/bin/bash
# =============================================================================
# 01 — liftOver mm39 anchor → mm10, then carve out a 200 kb locus FASTA.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

# 1. Build mm39 anchor + flank BED
ANCHOR_FLANK_MM39_BED="${TGT_REF_DIR}/anchor_flank_mm39.bed"
ANCHOR_FLANK_MM10_BED="${TGT_REF_DIR}/anchor_flank_mm10.bed"
ANCHOR_FLANK_UNMAP_BED="${TGT_REF_DIR}/anchor_flank_mm10.unmapped.bed"
ANCHOR_TIGHT_MM39_BED="${TGT_REF_DIR}/anchor_tight_mm39.bed"
ANCHOR_TIGHT_MM10_BED="${TGT_REF_DIR}/anchor_tight_mm10.bed"
ANCHOR_TIGHT_UNMAP_BED="${TGT_REF_DIR}/anchor_tight_mm10.unmapped.bed"

flank_start=$((ANCHOR_MM39_START - LOCUS_FLANK_BP))
flank_end=$((ANCHOR_MM39_END + LOCUS_FLANK_BP))
[ "${flank_start}" -lt 0 ] && flank_start=0

printf "%s\t%d\t%d\tsrrm3_locus_with_flanks\t.\t%s\n" \
  "${ANCHOR_MM39_CHROM}" "${flank_start}" "${flank_end}" "${ANCHOR_STRAND}" \
  > "${ANCHOR_FLANK_MM39_BED}"

printf "%s\t%d\t%d\tnovel_cassette_79bp\t.\t%s\n" \
  "${ANCHOR_MM39_CHROM}" "${ANCHOR_MM39_START}" "${ANCHOR_MM39_END}" "${ANCHOR_STRAND}" \
  > "${ANCHOR_TIGHT_MM39_BED}"

# 2. liftOver mm39 → mm10 (chain auto-fetched by the main pipeline already)
log_step "liftOver mm39 anchor + flanks → mm10"
liftOver "${ANCHOR_FLANK_MM39_BED}" "${CHAIN_MM39_TO_MM10}" \
         "${ANCHOR_FLANK_MM10_BED}" "${ANCHOR_FLANK_UNMAP_BED}"
liftOver "${ANCHOR_TIGHT_MM39_BED}" "${CHAIN_MM39_TO_MM10}" \
         "${ANCHOR_TIGHT_MM10_BED}" "${ANCHOR_TIGHT_UNMAP_BED}"

[ -s "${ANCHOR_FLANK_MM10_BED}" ] || { echo "FAIL: flank lift produced empty BED" >&2; exit 1; }
[ -s "${ANCHOR_TIGHT_MM10_BED}" ] || { echo "FAIL: tight anchor lift produced empty BED" >&2; exit 1; }

read -r LOCUS_CHR LOCUS_START LOCUS_END _LREST <"${ANCHOR_FLANK_MM10_BED}"
read -r ANCH_CHR ANCH_START ANCH_END _AREST <"${ANCHOR_TIGHT_MM10_BED}"

echo ""
echo "Locus (mm10): ${LOCUS_CHR}:${LOCUS_START}-${LOCUS_END}  ($((LOCUS_END - LOCUS_START)) bp)"
echo "Cassette (mm10): ${ANCH_CHR}:${ANCH_START}-${ANCH_END}  ($((ANCH_END - ANCH_START)) bp)"

# Save canonical anchor BED for downstream classification
cp -f "${ANCHOR_TIGHT_MM10_BED}" "${ANCHOR_MM10_BED}"

# 3. samtools faidx — extract 200 kb FASTA
log_step "samtools faidx mm10 → locus FASTA"
[ -f "${MM10_REF}.fai" ] || samtools faidx "${MM10_REF}"
samtools faidx "${MM10_REF}" \
  "${LOCUS_CHR}:$((LOCUS_START + 1))-${LOCUS_END}" \
  > "${LOCUS_FA}"
samtools faidx "${LOCUS_FA}"

# Save absolute locus origin for later coordinate translation in the
# CIGAR walker (the locus FASTA's coordinates are 1-based from start of
# the extracted region, not the original mm10 chromosome coordinates).
LOCUS_META="${TGT_REF_DIR}/locus_origin.tsv"
{
  printf "field\tvalue\n"
  printf "locus_fa_seqname\t%s:%d-%d\n" "${LOCUS_CHR}" "$((LOCUS_START + 1))" "${LOCUS_END}"
  printf "locus_chrom\t%s\n" "${LOCUS_CHR}"
  printf "locus_start_mm10\t%d\n" "${LOCUS_START}"
  printf "locus_end_mm10\t%d\n" "${LOCUS_END}"
  printf "anchor_start_mm10\t%d\n" "${ANCH_START}"
  printf "anchor_end_mm10\t%d\n" "${ANCH_END}"
  printf "anchor_strand\t%s\n" "${ANCHOR_STRAND}"
  # Position relative to locus FASTA (which is "0-based offset from extracted region start")
  printf "anchor_start_in_locus\t%d\n" "$((ANCH_START - LOCUS_START))"
  printf "anchor_end_in_locus\t%d\n" "$((ANCH_END - LOCUS_START))"
} > "${LOCUS_META}"

cat "${LOCUS_META}"

log_step "Done"
ls -lh "${TGT_REF_DIR}/"
