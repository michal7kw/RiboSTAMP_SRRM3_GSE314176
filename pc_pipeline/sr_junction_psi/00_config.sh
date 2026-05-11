# =============================================================================
# pc_pipeline/sr_junction_psi/00_config.sh
# Short-read junction-PSI pipeline. Tests the progenitor hypothesis (OPC,
# iDG-Immature) for the 79-bp Srrm3 cassette using rMATS-style logic on the
# matched short-read 10x library + cell-line bulk validation.
#
# Why this exists: the long-read targeted_psi/ pipeline measured ~0% PSI in
# 7 mature cell types, but couldn't test progenitors because OPC + iDG were
# filtered from the published long-read whitelist. The short-read library
# DOES contain those populations (175 + 221 = 396 cells), so we run a
# parallel junction-counting analysis here.
# =============================================================================

# Pull main config (BASE_DIR, MM10_REF, conda env, etc.)
THIS_DIR="$(dirname "$(readlink -f "$0")")"
MAIN_CONFIG="$(dirname "${THIS_DIR}")/00_config.sh"
# shellcheck disable=SC1090
source "${MAIN_CONFIG}"

# --- Override paths so we keep this analysis isolated -----------------------
export SR_DIR="${BASE_DIR}/sr_junction_psi"
export SR_BAM_DIR="${SR_DIR}/bams"
export SR_BAM_DIR_BULK="${SR_BAM_DIR}/cell_line_bulk"
export SR_BAM_DIR_10X="${SR_BAM_DIR}/short_read_10x"
export SR_PER_CELL_DIR="${SR_DIR}/per_cell"
export SR_RESULTS_DIR="${SR_DIR}/results"
export SR_LOG_DIR="${SR_DIR}/logs"

mkdir -p "${SR_BAM_DIR_BULK}" "${SR_BAM_DIR_10X}" \
         "${SR_PER_CELL_DIR}" "${SR_RESULTS_DIR}" "${SR_LOG_DIR}"

# --- Junction coordinates (lifted from mm39 anchor study to mm10) ----------
# From 90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md §3.4,
# applying the verified -28854 bp mm39→mm10 offset (see LIFTOVER_VERIFICATION).
#
# In mm10 (0-based, on chr5):
#   Inc1 junction: upstream-exon-end → cassette-start  (135869294 → 135869720)
#   Inc2 junction: cassette-end → downstream-exon-start (135869798 → 135873080)
#   Skip junction: upstream-exon-end → downstream-exon-start (135869294 → 135873080)
#
# These are the three splice-junction coordinates we count reads at.
export SR_CHROM="chr5"
export SR_INC1_LEFT=135869294    # upstream exon end on +
export SR_INC1_RIGHT=135869720   # cassette start on +
export SR_INC2_LEFT=135869798    # cassette end on +
export SR_INC2_RIGHT=135873080   # downstream exon start on +
export SR_SKIP_LEFT="${SR_INC1_LEFT}"
export SR_SKIP_RIGHT="${SR_INC2_RIGHT}"

# Tolerance window for a CIGAR N op to "match" each junction. rMATS uses
# ±0 (exact); we allow ±5 bp to absorb minor alignment slop. A read counts
# as supporting that junction if its CIGAR has an N op whose left edge is
# within ±SR_JUNC_TOL of the junction's left coordinate AND whose right edge
# is within ±SR_JUNC_TOL of the right coordinate.
export SR_JUNC_TOL=5

# --- Sample accessions (resolved 2026-04-28 via pysradb) -------------------
#
# Short-read 10x scRNA-seq (matched to long-read library — same 3 mice).
# These contain OPC + iDG-Immature populations (the ones missing from the
# long-read whitelist).
#
# IMPORTANT: each GSM has 2 SRR runs (paired-end Element AVITI). The data is
# RAW FASTQ — to get CB-tagged BAMs we need to align with cellranger or
# starsolo first. ~150 GB raw FASTQ per mouse × 3 = ~450-500 GB total.
# This is a multi-day pipeline, not a single-script run.
export SR_10X_GSMS=("GSM9380796" "GSM9380797" "GSM9380798")
declare -A SR_10X_GSM_TO_SRRS=(
    [GSM9380796]="SRR36480583 SRR36480584"   # mouse1 (rep 1) — paired Element AVITI runs
    [GSM9380797]="SRR36480581 SRR36480582"   # mouse2 (rep 2)
    [GSM9380798]="SRR36480579 SRR36480580"   # mouse3 (rep 3)
)
# Note: SRA library_strategy=RNA-Seq, library_source=TRANSCRIPTOMIC SINGLE CELL,
# instrument=Element AVITI. The raw FASTQ needs to be processed through
# cellranger / starsolo to produce a CB-tagged BAM before our junction-PSI
# calculator can do per-cell aggregation. See README.md §"Week 1 — 10x reality".

# Cell-line Ribo-STAMP / bulk RNA-seq (validation set — should give ~57% PSI
# in NT samples, matching the lab's bulk anchor study).
# These are Illumina NovaSeq 6000, single-end, ~2 GB SRA each (~18 GB total).
# Tractable in 1-2 hours of compute including download + STAR alignment.
export SR_BULK_GSMS=("GSM8465528" "GSM8465529" "GSM8465530"
                     "GSM8465531" "GSM8465532" "GSM8465533"
                     "GSM8465534" "GSM8465535" "GSM8465536")
declare -A SR_BULK_GSM_TO_SRRS=(
    # NT samples (rep 1-3) — WT condition, expected PSI ≈ 57%
    [GSM8465528]="SRR30278197 SRR30278198"   # NTrep1
    [GSM8465529]="SRR30278195 SRR30278196"   # NTrep2
    [GSM8465530]="SRR30278193 SRR30278194"   # NTrep3
    # B15 samples (rep 1-3) — 15 min BDNF stimulation
    [GSM8465531]="SRR30278191 SRR30278192"   # B15rep1
    [GSM8465532]="SRR30278189 SRR30278190"   # B15rep2
    [GSM8465533]="SRR30278187 SRR30278188"   # B15rep3
    # B60 samples (rep 1-3) — 60 min BDNF stimulation
    [GSM8465534]="SRR30278185 SRR30278186"   # B60rep1
    [GSM8465535]="SRR30278183 SRR30278184"   # B60rep2
    [GSM8465536]="SRR30278181 SRR30278182"   # B60rep3
)
# Note: cell-line BAMs need STAR alignment before junction-PSI counting.
# The bulk samples are tractable: ~30 min download + ~1 hr STAR + ~5 min PSI.

# --- Per-cell aggregation parameters ----------------------------------------
# Cells with very few reads at the locus can't vote meaningfully. Filter to
# cells with at least N reads at the cassette region.
export SR_MIN_READS_PER_CELL=3
export SR_MIN_READS_PER_CLUSTER=50

# CB extraction (short-read 10x uses CB tag in the BAM, simpler than long-read)
export SR_CB_TAG="CB"   # Cell barcode tag
export SR_UB_TAG="UB"   # UMI tag

# --- Thread budget ----------------------------------------------------------
export SR_NUM_THREADS="${SR_NUM_THREADS:-16}"

# --- Helpers ---------------------------------------------------------------
log_step() {
  printf "\n=== [%s] %s ===\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}
