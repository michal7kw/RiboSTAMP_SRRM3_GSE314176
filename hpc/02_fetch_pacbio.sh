#!/bin/bash
#SBATCH --job-name=ribo_02_pacbio
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --array=0-2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/02_pacbio_%A_%a.err
#SBATCH --output=./logs/02_pacbio_%A_%a.out

# =============================================================================
# Download three PacBio Revio HiFi runs for the Hippocampus dataset
# (the only PACBIO_SMRT platform runs under BioProject PRJNA1389908).
#
# Sample metadata (from SRA Run Selector, as of 2026-03):
#   SRR36480452  GSM9380801  115 GB / 35.7 Gbases  batch 20241229
#   SRR36480453  GSM9380800  127 GB / 40.2 Gbases  batch 20241229
#   SRR36480454  GSM9380799  131 GB / 40.4 Gbases  batch 20241228
# All three: dHC injection with 800 nL EF1a-Ribo STAMP virus, C57BL/6J,
# MAS-seq concatemers (~13,500 bp avg spot length — need skera split next).
#
# Total raw download: ~374 GB compressed, more after sam-dump.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
RAW_DIR="${BASE_DIR}/data/raw_pacbio"
mkdir -p "${RAW_DIR}" "$(pwd)/logs"

HPC_SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
ACC_LIST="${HPC_SCRIPT_DIR}/pacbio_SRR_Acc_List.txt"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

# Select this task's accession (one SRR per array task)
mapfile -t SRRS < <(grep -v '^\s*$' "${ACC_LIST}")
SRR="${SRRS[$SLURM_ARRAY_TASK_ID]}"

echo "============================================"
echo "prefetch + sam-dump  — ${SRR}"
echo "Start: $(date)"
echo "============================================"

OUT_BAM="${RAW_DIR}/${SRR}.unaligned.bam"
if [ -s "${OUT_BAM}" ]; then
  echo "${OUT_BAM} already exists and is non-empty; skipping."
  exit 0
fi

# 1) prefetch to a per-run directory (sra-tools creates ${SRR}/${SRR}.sra inside).
#    --max-size increased because PacBio runs can be >100 GB.
prefetch \
    --output-directory "${RAW_DIR}" \
    --max-size u \
    --resume yes \
    "${SRR}"

SRA_PATH="${RAW_DIR}/${SRR}/${SRR}.sra"
if [ ! -f "${SRA_PATH}" ]; then
  # Some sra-tools versions flatten the layout. Try the direct path.
  SRA_PATH="${RAW_DIR}/${SRR}.sra"
fi
if [ ! -f "${SRA_PATH}" ]; then
  echo "ERROR: prefetch succeeded but ${SRR}.sra is not where expected." >&2
  find "${RAW_DIR}" -name "${SRR}*.sra" >&2
  exit 1
fi

# 2) sam-dump -> BAM. We use sam-dump (not fasterq-dump) because PacBio HiFi
#    reads carry CCS/kinetic BAM tags that FASTQ would drop; skera + pbmm2 need them.
echo "sam-dump ${SRR} -> ${OUT_BAM}"
sam-dump --unaligned "${SRA_PATH}" \
  | samtools view -b - \
  > "${OUT_BAM}"

# 3) sanity check: the file should be > 1 GB for any PacBio HiFi run.
BYTES=$(stat -c%s "${OUT_BAM}")
if [ "${BYTES}" -lt 1000000000 ]; then
  echo "WARNING: ${OUT_BAM} is only ${BYTES} bytes — likely a partial dump. Re-run." >&2
fi

echo "=== Done: ${SRR} ($(numfmt --to=iec "${BYTES}"))  $(date) ==="
