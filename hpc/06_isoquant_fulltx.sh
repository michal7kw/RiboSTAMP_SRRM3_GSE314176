#!/bin/bash
#SBATCH --job-name=ribo_06_isoq_full
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/06_isoq_full_%j.err
#SBATCH --output=./logs/06_isoq_full_%j.out

# =============================================================================
# ON-DEMAND full-transcriptome IsoQuant run.
#
# Do NOT submit unless:
#   (a) hpc/05_isoquant_targeted.sh flagged a novel isoform, AND
#   (b) you want to confirm locus-subsetting didn't bias read support.
#
# Verification criterion (see plan file): per-cell novel-isoform read count
# Spearman r > 0.95 between this run and the targeted run. Below that ->
# targeted subsetting introduced bias; reject the targeted result.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
ALIGN_DIR="${BASE_DIR}/data/aligned_mm10"
OUT_DIR="${BASE_DIR}/data/isoquant_fulltx"
mkdir -p "${OUT_DIR}" "$(pwd)/logs"

MM10_REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/mm10.fa"
GENCODE_GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/gencode.vM25.annotation.gtf"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

isoquant \
    --reference "${MM10_REF}" \
    --genedb "${GENCODE_GTF}" \
    --bam "${ALIGN_DIR}"/*.mm10.bam \
    --data_type pacbio_ccs \
    --output "${OUT_DIR}" \
    --threads "${SLURM_NTASKS}" \
    --sqanti_output \
    --count_exons \
    --read_group tag:CB

echo "=== Done: $(date) ==="
