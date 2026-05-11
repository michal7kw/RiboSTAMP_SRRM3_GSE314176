# =============================================================================
# pc_pipeline/targeted_psi/00_config.sh
# Sourced by every script in this folder. Inherits BASE_DIR & friends from
# the main pipeline's 00_config.sh, then overrides paths so we write to a
# completely separate data/targeted_psi/ tree.
# =============================================================================

# Pull the main config first (BASE_DIR, MM10_REF, GENCODE_GTF, conda env, etc).
THIS_DIR="$(dirname "$(readlink -f "$0")")"
MAIN_CONFIG="$(dirname "${THIS_DIR}")/00_config.sh"
# shellcheck disable=SC1090
source "${MAIN_CONFIG}"

# --- Override paths so we don't touch the main pipeline's outputs ----------
export TGT_DIR="${BASE_DIR}/targeted_psi"
export TGT_REF_DIR="${TGT_DIR}/reference"
export TGT_ALIGN_DIR="${TGT_DIR}/aligned_locus"
export TGT_PERREAD_DIR="${TGT_DIR}/per_read"
export TGT_RESULTS_DIR="${TGT_DIR}/results"
export TGT_LOG_DIR="${TGT_DIR}/logs"

mkdir -p "${TGT_REF_DIR}" "${TGT_ALIGN_DIR}" "${TGT_PERREAD_DIR}" \
         "${TGT_RESULTS_DIR}" "${TGT_LOG_DIR}"

# Locus FASTA + anchor BEDs (mm10) — produced by 01_extract_locus.sh
export LOCUS_FA="${TGT_REF_DIR}/srrm3_locus_mm10.fa"
export ANCHOR_MM10_BED="${TGT_REF_DIR}/anchor_mm10.bed"

# --- The lab short-read anchor in mm39 ---------------------------------------
# From 90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md
# 79-bp cassette, negative strand
export ANCHOR_MM39_CHROM="chr5"
export ANCHOR_MM39_START="135898574"
export ANCHOR_MM39_END="135898652"
export ANCHOR_LENGTH=79
export ANCHOR_STRAND="-"

# Context window around the cassette (in mm39 coords) for the locus FASTA.
# Originally 200 kb — empirically too wide: caught ~3M spurious-seed reads
# per BAM in adjacent genes + repeats, which slowed Phase B (pbmm2 splice)
# from estimated 5-15 min to actual 4.5 hours. Shrinking to 30 kb covers
# SRRM3's longest known introns (~25 kb) with a small safety margin and
# excludes flanking genes. See LEARNING_07 for the postmortem.
export LOCUS_FLANK_BP="${LOCUS_FLANK_BP:-30000}"

# --- Thread budget ----------------------------------------------------------
# Empirically (after a 2h failed run with 3-parallel × 5 threads) this
# workload is BOTH CPU- AND memory-bandwidth-sensitive AND read-I/O on the
# WSL2 9P filesystem suffers under concurrent streams. Sequential 1 SRR at
# a time × 16 threads gets full memory bandwidth + single-stream sequential
# read on /mnt/e/, which is much faster wall-clock than 3-in-parallel.
export TGT_NUM_THREADS="${TGT_NUM_THREADS:-16}"
export TGT_PARALLEL_SRRS="${TGT_PARALLEL_SRRS:-1}"

# --- Barcode extraction parameters ------------------------------------------
# 10x Chromium Next GEM 3' v3.4 (per the SRA metadata) uses a 16-bp cell
# barcode + 12-bp UMI. The CB sits adjacent to the polyA.
export CB_LEN=16
export UMI_LEN=12
export POLYA_MIN=8         # minimum run of A/T to anchor barcode position
export BARCODE_MAX_HAMMING=1   # max mismatches allowed when matching to whitelist

# --- Helpers ---------------------------------------------------------------
log_step() {
  printf "\n=== [%s] %s ===\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}
