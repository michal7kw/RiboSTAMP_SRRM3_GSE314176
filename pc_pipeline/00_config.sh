# =============================================================================
# pc_pipeline/00_config.sh — central config sourced by every pc_pipeline script.
#
# Edit the values below to match your machine, then every other script in this
# folder uses these settings. This is the only file that should need changing.
# =============================================================================

# --- Where data lives -------------------------------------------------------
# Default: E: drive in Windows ↔ /mnt/e/ in WSL.
# Trade-off note: /mnt/e/ goes through WSL's 9P translation layer and is
# ~2-3× slower than WSL-native ext4. For typical bioinformatics I/O this is
# fine. If you want max speed, point BASE_DIR at ~/ribodata/ (inside WSL home,
# invisible from Windows Explorer but on native ext4).
export BASE_DIR="${BASE_DIR:-/mnt/e/RiboSTAMP_SRRM3_GSE314176/data}"

# --- Conda env --------------------------------------------------------------
# Whatever 'conda info --base' returns on your WSL install. Common locations:
#   ~/miniconda3
#   ~/anaconda3
#   /opt/miniconda3
export CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"

# Name of the env (created from ../env/longread_iso.yml during setup).
export LONGREAD_ENV="${LONGREAD_ENV:-longread_iso}"

# SQANTI3 lives in its own env (it is not in bioconda — see setup_sqanti3.sh).
export SQANTI3_ENV="${SQANTI3_ENV:-sqanti3}"
export SQANTI3_HOME="${SQANTI3_HOME:-$HOME/SQANTI3}"

# --- Reference paths --------------------------------------------------------
# Same logical layout as the HPC version, just rooted under BASE_DIR.
export REF_DIR="${REF_DIR:-${BASE_DIR}/reference/mm10}"
export MM10_REF="${MM10_REF:-${REF_DIR}/mm10.fa}"
export GENCODE_GTF="${GENCODE_GTF:-${REF_DIR}/gencode.vM25.annotation.gtf}"
export CHAIN_DIR="${CHAIN_DIR:-${BASE_DIR}/reference/chains}"
export CHAIN_MM10_TO_MM39="${CHAIN_MM10_TO_MM39:-${CHAIN_DIR}/mm10ToMm39.over.chain.gz}"
export CHAIN_MM39_TO_MM10="${CHAIN_MM39_TO_MM10:-${CHAIN_DIR}/mm39ToMm10.over.chain.gz}"

# --- Compute knobs ----------------------------------------------------------
# i5-14600K = 14 cores / 20 threads. We default to 16 threads per tool and
# keep 4 free for the OS / I/O. Adjust for your machine.
export NUM_THREADS="${NUM_THREADS:-16}"

# --- Which SRR runs to process ----------------------------------------------
# Default: all three Revio runs. For quick_test.sh, override to just the
# smallest (SRR36480452, ~115 GB) to validate the pipeline end-to-end.
export SRR_LIST="${SRR_LIST:-SRR36480452 SRR36480453 SRR36480454}"

# --- Internal: derived paths used by the scripts ----------------------------
# (Mirror the layout of the HPC version so docs/LEARNING_*.md still apply.)
export PROC_DIR="${BASE_DIR}/processed"
export META_DIR="${BASE_DIR}/metadata"
export RAW_DIR="${BASE_DIR}/raw_pacbio"
export SPLIT_DIR="${BASE_DIR}/skera_split"
export ALIGN_DIR="${BASE_DIR}/aligned_mm10"
export ANCHOR_DIR="${BASE_DIR}/anchor_coverage"
export ISOQ_TGT_DIR="${BASE_DIR}/isoquant_targeted"
export ISOQ_FULL_DIR="${BASE_DIR}/isoquant_fulltx"
export SQANTI_DIR="${BASE_DIR}/sqanti3"
export LIFT_DIR="${BASE_DIR}/liftover"
export SUBSET_DIR="${BASE_DIR}/srrm3_subset"
export PC_LOG_DIR="${BASE_DIR}/logs"

# Path to the Srrm3 locus BED (lives with the HPC scripts; reused as-is)
export LOCUS_BED="$(dirname "$(readlink -f "$0")")/../hpc/srrm3_locus.bed"

# --- Helpers ----------------------------------------------------------------
mkdir -p "${PROC_DIR}" "${META_DIR}" "${RAW_DIR}" "${SPLIT_DIR}" \
         "${ALIGN_DIR}" "${ANCHOR_DIR}" "${ISOQ_TGT_DIR}" \
         "${ISOQ_FULL_DIR}" "${SQANTI_DIR}" "${LIFT_DIR}" \
         "${SUBSET_DIR}" "${PC_LOG_DIR}"

activate_env() {
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "${LONGREAD_ENV}"
}

log_step() {
  printf "\n=== [%s] %s ===\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}
