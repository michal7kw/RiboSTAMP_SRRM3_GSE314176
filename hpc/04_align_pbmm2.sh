#!/bin/bash
#SBATCH --job-name=ribo_04_pbmm2
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/04_pbmm2_%A_%a.err
#SBATCH --output=./logs/04_pbmm2_%A_%a.out

# =============================================================================
# pbmm2 alignment of post-skera S-reads to mm10 + Gencode vM25.
#
# Build choice: mm10 matches the authors' alignments so their cell-barcode-to-
# read correspondences (obs.csv) remain valid. mm39 cross-reference happens
# later in 08_liftover_mm10_to_mm39.sh.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
SPLIT_DIR="${BASE_DIR}/data/skera_split"
ALIGN_DIR="${BASE_DIR}/data/aligned_mm10"
mkdir -p "${ALIGN_DIR}" "$(pwd)/logs"

HPC_SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
ACC_LIST="${HPC_SCRIPT_DIR}/pacbio_SRR_Acc_List.txt"

# TODO: provision mm10 reference on /beegfs if not already present.
# Download: rsync -avz rsync://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz .
MM10_REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/mm10.fa"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

mapfile -t SRRS < <(grep -v '^\s*$' "${ACC_LIST}")
SRR="${SRRS[$SLURM_ARRAY_TASK_ID]}"

IN_BAM="${SPLIT_DIR}/${SRR}.segmented.bam"
OUT_BAM="${ALIGN_DIR}/${SRR}.mm10.bam"

echo "============================================"
echo "pbmm2 align — ${SRR}"
echo "Start: $(date)"
echo "============================================"

if [ ! -f "${IN_BAM}" ]; then
  echo "ERROR: missing ${IN_BAM} — did 03_skera_split.sh complete for ${SRR}?" >&2
  exit 1
fi

if [ -s "${OUT_BAM}" ]; then
  echo "Already aligned; skipping."
  exit 0
fi

pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --num-threads "${SLURM_NTASKS}" \
    "${MM10_REF}" \
    "${IN_BAM}" \
    "${OUT_BAM}"

samtools index "${OUT_BAM}"

echo "=== Done: ${SRR}  $(date) ==="
