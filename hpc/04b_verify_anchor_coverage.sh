#!/bin/bash
#SBATCH --job-name=ribo_04b_anchor
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/04b_anchor_%j.err
#SBATCH --output=./logs/04b_anchor_%j.out

# =============================================================================
# Tool-independent positive-control check for the novel Srrm3 cassette.
#
# Rationale (see docs/LEARNING_06_novelty_layered_defense.md):
#   The lab's short-read pipeline identified a 79-bp novel cassette at
#   chr5:135,898,574-135,898,652 (mm39). Before trusting any IsoQuant output,
#   we want to verify *at the BAM level* that:
#     (a) the alignment actually covers the cassette coordinates, and
#     (b) some reads exhibit the splice-junction pattern that defines the
#         novel cassette (skipping into and out of the 79-bp region).
#
# Why "tool-independent": if IsoQuant later produces no novel transcript for
# Srrm3, we need to know whether the failure is in alignment, in splice
# detection, or in IsoQuant's discovery threshold — three very different
# debugging paths. samtools + CIGAR parsing answers (a) and (b) with no
# dependency on IsoQuant's internal logic.
#
# Run order: after 04_align_pbmm2.sh, BEFORE 05_isoquant_targeted.sh.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
ALIGN_DIR="${BASE_DIR}/data/aligned_mm10"
OUT_DIR="${BASE_DIR}/data/anchor_coverage"
mkdir -p "${OUT_DIR}" "$(pwd)/logs"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

# --- Step 1: lift the mm39 anchor to mm10 -----------------------------------
# We aligned to mm10 (matches authors). The anchor was identified in mm39.
# liftOver here is the mirror of step 08 (which goes the other direction).

CHAIN="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm39ToMm10.over.chain.gz"
if [ ! -f "${CHAIN}" ]; then
  echo "Fetching mm39ToMm10 chain (one-time)..."
  mkdir -p "$(dirname "${CHAIN}")"
  wget -O "${CHAIN}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm10.over.chain.gz"
fi

ANCHOR_MM39="${OUT_DIR}/anchor_mm39.bed"
ANCHOR_MM10="${OUT_DIR}/anchor_mm10.bed"
ANCHOR_UNMAPPED="${OUT_DIR}/anchor_mm10.unmapped.bed"

# Anchor BED with 1 kb flanks for context-tolerant coverage check
printf "chr5\t135897574\t135899652\tnovel_srrm3_cassette_with_1kb_flanks\t.\t-\n" \
  > "${ANCHOR_MM39}"

liftOver "${ANCHOR_MM39}" "${CHAIN}" "${ANCHOR_MM10}" "${ANCHOR_UNMAPPED}"

if [ ! -s "${ANCHOR_MM10}" ]; then
  echo "ERROR: liftOver failed to map the anchor from mm39 to mm10." >&2
  echo "Unmapped reasons:" >&2
  cat "${ANCHOR_UNMAPPED}" >&2
  exit 2
fi

read -r CHR ANCHOR_START ANCHOR_END _REST <"${ANCHOR_MM10}"
ANCHOR_REGION="${CHR}:${ANCHOR_START}-${ANCHOR_END}"
echo "Anchor lifted to mm10: ${ANCHOR_REGION}"

# Tighter coordinates for the cassette ITSELF (no flanks) — used for the
# junction-pattern check below. The 79-bp cassette in mm10 may be a few bp
# shifted from mm39; recompute by lifting the strict anchor too.
STRICT_MM39="${OUT_DIR}/anchor_strict_mm39.bed"
STRICT_MM10="${OUT_DIR}/anchor_strict_mm10.bed"
printf "chr5\t135898574\t135898652\tnovel_cassette_strict\t.\t-\n" > "${STRICT_MM39}"
liftOver "${STRICT_MM39}" "${CHAIN}" "${STRICT_MM10}" "${OUT_DIR}/anchor_strict_mm10.unmapped.bed"

read -r _ STRICT_START STRICT_END _ <"${STRICT_MM10}"
echo "Strict cassette in mm10: ${CHR}:${STRICT_START}-${STRICT_END} (length $((STRICT_END - STRICT_START)) bp)"

# --- Step 2: per-BAM coverage and junction checks ----------------------------

REPORT="${OUT_DIR}/anchor_coverage_report.tsv"
{
  printf "sample\tregion_with_flanks\treads_in_region\tspliced_reads_in_region\treads_with_junction_pattern\tverdict\n"
} > "${REPORT}"

ANY_FAIL=0
for BAM in "${ALIGN_DIR}"/*.mm10.bam; do
  [ -f "${BAM}" ] || continue
  NAME="$(basename "${BAM}" .mm10.bam)"
  echo "--- ${NAME} ---"

  # (a) total reads spanning the lifted anchor + flanks
  N_READS=$(samtools view -c "${BAM}" "${ANCHOR_REGION}")

  # (b) spliced reads only (CIGAR contains 'N' = ref-skip, i.e. exon-exon junction)
  N_SPLICED=$(samtools view "${BAM}" "${ANCHOR_REGION}" | awk '$6 ~ /N/' | wc -l)

  # (c) reads whose CIGAR contains a splice junction (N segment) within
  # +/- 50 bp of the strict cassette boundaries — i.e. the read appears to
  # use a splice site that DEFINES the novel cassette. This requires parsing
  # CIGAR + computing reference positions.
  N_JUNC=$(samtools view "${BAM}" "${ANCHOR_REGION}" \
    | awk -v cs="${STRICT_START}" -v ce="${STRICT_END}" -v tol=50 '
      {
        cigar = $6; pos = $4   # 1-based leftmost ref position
        # Walk the CIGAR; for each N op record the (ref_start, ref_end) of the skip
        ref = pos
        n = length(cigar)
        op_start = 1
        for (i = 1; i <= n; i++) {
          c = substr(cigar, i, 1)
          if (c ~ /[A-Z=]/) {
            op = c
            len = substr(cigar, op_start, i - op_start) + 0
            op_start = i + 1
            if (op == "N") {
              skip_start = ref
              skip_end   = ref + len - 1
              # Hits a splice site at the upstream cassette boundary?
              if ( (skip_end+1 >= cs - tol && skip_end+1 <= cs + tol) ||
                   (skip_start  >= ce - tol && skip_start  <= ce + tol) ) {
                hit = 1
                break
              }
              ref += len
            } else if (op == "M" || op == "D" || op == "=" || op == "X") {
              ref += len
            }
            # I, S, H, P do not consume reference
          }
        }
        if (hit) print
        hit = 0
      }' | wc -l)

  # Verdict: at minimum we want >=10 reads in region AND >=1 with the cassette
  # junction pattern. Tunable; see LEARNING_06.
  if [ "${N_READS}" -lt 10 ]; then
    VERDICT="FAIL_LOW_COVERAGE"
    ANY_FAIL=1
  elif [ "${N_JUNC}" -lt 1 ]; then
    VERDICT="WARN_NO_CASSETTE_JUNCTION"
  else
    VERDICT="PASS"
  fi

  printf "%s\t%s\t%d\t%d\t%d\t%s\n" \
    "${NAME}" "${ANCHOR_REGION}" "${N_READS}" "${N_SPLICED}" "${N_JUNC}" "${VERDICT}" \
    >> "${REPORT}"

  echo "  total reads:                ${N_READS}"
  echo "  spliced reads:              ${N_SPLICED}"
  echo "  reads with cassette splice: ${N_JUNC}"
  echo "  verdict:                    ${VERDICT}"
done

echo ""
echo "=== Report: ${REPORT} ==="
column -t -s $'\t' "${REPORT}"

if [ "${ANY_FAIL}" -ne 0 ]; then
  echo ""
  echo "FAIL: at least one sample has < 10 reads at the anchor locus." >&2
  echo "This means alignment didn't cover the Srrm3 cassette region. Either" >&2
  echo "the mm10 reference is wrong or skera/pbmm2 had a problem. Do NOT" >&2
  echo "proceed to 05_isoquant_targeted.sh until this is resolved." >&2
  exit 1
fi

echo ""
echo "All samples have coverage at the anchor locus. Proceed to 05_isoquant_targeted.sh."
echo "If verdict is WARN_NO_CASSETTE_JUNCTION on all samples, that is a real biological"
echo "finding (no PacBio reads support the novel cassette in this dataset). See"
echo "docs/LEARNING_06_novelty_layered_defense.md for interpretation."
echo "=== Done: $(date) ==="
