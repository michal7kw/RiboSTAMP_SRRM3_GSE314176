#!/bin/bash
# =============================================================================
# pc_pipeline/04b — BLOCKING tool-independent positive control.
# Verify the alignment actually covers the cassette + uses its splice junctions.
# See docs/LEARNING_06_novelty_layered_defense.md.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

mkdir -p "${ANCHOR_DIR}"

# --- Step 1: lift mm39 anchor to mm10 ---------------------------------------
if [ ! -f "${CHAIN_MM39_TO_MM10}" ]; then
  log_step "Fetching mm39ToMm10 chain (one-time)"
  mkdir -p "${CHAIN_DIR}"
  wget -O "${CHAIN_MM39_TO_MM10}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm10.over.chain.gz"
fi

ANCHOR_MM39="${ANCHOR_DIR}/anchor_mm39.bed"
ANCHOR_MM10="${ANCHOR_DIR}/anchor_mm10.bed"
printf "chr5\t135897574\t135899652\tnovel_srrm3_cassette_with_1kb_flanks\t.\t-\n" > "${ANCHOR_MM39}"
liftOver "${ANCHOR_MM39}" "${CHAIN_MM39_TO_MM10}" "${ANCHOR_MM10}" "${ANCHOR_DIR}/anchor_mm10.unmapped.bed"

if [ ! -s "${ANCHOR_MM10}" ]; then
  echo "ERROR: liftOver failed for the anchor." >&2
  cat "${ANCHOR_DIR}/anchor_mm10.unmapped.bed" >&2
  exit 2
fi

read -r CHR ANCHOR_START ANCHOR_END _ <"${ANCHOR_MM10}"
ANCHOR_REGION="${CHR}:${ANCHOR_START}-${ANCHOR_END}"
echo "Anchor lifted to mm10: ${ANCHOR_REGION}"

STRICT_MM39="${ANCHOR_DIR}/anchor_strict_mm39.bed"
STRICT_MM10="${ANCHOR_DIR}/anchor_strict_mm10.bed"
printf "chr5\t135898574\t135898652\tnovel_cassette_strict\t.\t-\n" > "${STRICT_MM39}"
liftOver "${STRICT_MM39}" "${CHAIN_MM39_TO_MM10}" "${STRICT_MM10}" "${ANCHOR_DIR}/anchor_strict_mm10.unmapped.bed"
read -r _ STRICT_START STRICT_END _ <"${STRICT_MM10}"
echo "Strict cassette in mm10: ${CHR}:${STRICT_START}-${STRICT_END}"

# --- Step 2: per-BAM coverage + junction-pattern checks ----------------------
REPORT="${ANCHOR_DIR}/anchor_coverage_report.tsv"
printf "sample\tregion\treads_in_region\tspliced_reads\treads_with_cassette_junction\tverdict\n" > "${REPORT}"

ANY_FAIL=0
for BAM in "${ALIGN_DIR}"/*.mm10.bam; do
  [ -f "${BAM}" ] || continue
  NAME="$(basename "${BAM}" .mm10.bam)"
  log_step "${NAME}: anchor coverage check"

  N_READS=$(samtools view -c "${BAM}" "${ANCHOR_REGION}")
  N_SPLICED=$(samtools view "${BAM}" "${ANCHOR_REGION}" | awk '$6 ~ /N/' | wc -l)
  N_JUNC=$(samtools view "${BAM}" "${ANCHOR_REGION}" \
    | awk -v cs="${STRICT_START}" -v ce="${STRICT_END}" -v tol=50 '
      {
        cigar = $6; ref = $4; n = length(cigar); op_start = 1; hit = 0
        for (i = 1; i <= n; i++) {
          c = substr(cigar, i, 1)
          if (c ~ /[A-Z=]/) {
            op = c
            len = substr(cigar, op_start, i - op_start) + 0
            op_start = i + 1
            if (op == "N") {
              skip_start = ref; skip_end = ref + len - 1
              if ( (skip_end+1 >= cs - tol && skip_end+1 <= cs + tol) ||
                   (skip_start  >= ce - tol && skip_start  <= ce + tol) ) hit = 1
              ref += len
            } else if (op == "M" || op == "D" || op == "=" || op == "X") {
              ref += len
            }
          }
        }
        if (hit) print
      }' | wc -l)

  if [ "${N_READS}" -lt 10 ]; then
    VERDICT="FAIL_LOW_COVERAGE"; ANY_FAIL=1
  elif [ "${N_JUNC}" -lt 1 ]; then
    VERDICT="WARN_NO_CASSETTE_JUNCTION"
  else
    VERDICT="PASS"
  fi

  printf "%s\t%s\t%d\t%d\t%d\t%s\n" \
    "${NAME}" "${ANCHOR_REGION}" "${N_READS}" "${N_SPLICED}" "${N_JUNC}" "${VERDICT}" \
    >> "${REPORT}"
  echo "  reads=${N_READS}  spliced=${N_SPLICED}  cassette-junction=${N_JUNC}  -> ${VERDICT}"
done

echo ""
column -t -s $'\t' "${REPORT}"

if [ "${ANY_FAIL}" -ne 0 ]; then
  echo ""
  echo "FAIL: at least one sample has < 10 reads at the anchor locus." >&2
  echo "Do NOT proceed to step 05 until alignment coverage is fixed." >&2
  exit 1
fi

log_step "Anchor coverage check passed for all samples."
