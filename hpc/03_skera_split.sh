#!/bin/bash
#SBATCH --job-name=ribo_03_skera
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/03_skera_%A_%a.err
#SBATCH --output=./logs/03_skera_%A_%a.out

# =============================================================================
# Skera split — deconcatenate MAS-seq arrays into individual S-reads.
#
# The Revio runs (~13,500 bp average spot length) are concatemers: multiple
# cDNAs joined end-to-end by the MAS adapter set ("masSeq"). skera reads the
# adapter positions from the PacBio BAM and splits them into separate reads,
# one per original cDNA, each carrying the source concatemer's cell barcode.
#
# Output BAM is what pbmm2 actually aligns. Without this step, pbmm2 would
# try to align 13.5 kb concatemers as single transcripts — garbage.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
RAW_DIR="${BASE_DIR}/data/raw_pacbio"
SPLIT_DIR="${BASE_DIR}/data/skera_split"
mkdir -p "${SPLIT_DIR}" "$(pwd)/logs"

HPC_SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
ACC_LIST="${HPC_SCRIPT_DIR}/pacbio_SRR_Acc_List.txt"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

mapfile -t SRRS < <(grep -v '^\s*$' "${ACC_LIST}")
SRR="${SRRS[$SLURM_ARRAY_TASK_ID]}"

IN_BAM="${RAW_DIR}/${SRR}.unaligned.bam"
OUT_BAM="${SPLIT_DIR}/${SRR}.segmented.bam"

echo "============================================"
echo "skera split — ${SRR}"
echo "Start:  $(date)"
echo "Input:  ${IN_BAM}"
echo "Output: ${OUT_BAM}"
echo "============================================"

if [ ! -f "${IN_BAM}" ]; then
  echo "ERROR: missing ${IN_BAM} — did 02_fetch_pacbio.sh complete for ${SRR}?" >&2
  exit 1
fi

if [ -s "${OUT_BAM}" ]; then
  echo "Already split; skipping."
  exit 0
fi

# Adapter set: MAS-Iso-seq v1 uses mas16_primers.fasta (shipped with skera).
# The library prep kit is confirmed as Kinnex full-length RNA (aka MAS-seq for 10x)
# from the paper's description (Ribo-STAMP + MAS-ISO-seq). If skera complains
# about adapter mismatches, check the exact kit version in the paper's methods
# and swap in mas8_primers.fasta or the kit-specific adapter FASTA.
ADAPTER_FASTA="$(dirname "$(which skera)")/../share/skera/mas16_primers.fasta"
if [ ! -f "${ADAPTER_FASTA}" ]; then
  # Fall back to a common conda location layout.
  ADAPTER_FASTA="$(python -c 'import skera, os; print(os.path.join(os.path.dirname(skera.__file__), "mas16_primers.fasta"))' 2>/dev/null || true)"
fi
if [ ! -f "${ADAPTER_FASTA}" ]; then
  echo "ERROR: could not locate MAS adapter FASTA. Run 'skera --help' and set ADAPTER_FASTA manually." >&2
  exit 1
fi
echo "Adapter FASTA: ${ADAPTER_FASTA}"

skera split \
    --num-threads "${SLURM_NTASKS}" \
    "${IN_BAM}" \
    "${ADAPTER_FASTA}" \
    "${OUT_BAM}"

# skera writes a summary.json next to the output — log read counts so we can
# sanity-check that split reads are ~N x concatemer reads (typical MAS-16 = 16x).
SUMMARY="${OUT_BAM%.bam}.summary.json"
if [ -f "${SUMMARY}" ]; then
  echo "--- skera summary ---"
  cat "${SUMMARY}"
fi

echo "=== Done: ${SRR}  $(date) ==="
