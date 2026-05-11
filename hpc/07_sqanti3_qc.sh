#!/bin/bash
#SBATCH --job-name=ribo_07_sqanti
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/07_sqanti_%j.err
#SBATCH --output=./logs/07_sqanti_%j.out

# =============================================================================
# Classify IsoQuant output with SQANTI3.
#
# Produces FSM / ISM / NIC / NNC / genic / intergenic / antisense / fusion
# categories + flags RT-switching and intra-priming artifacts. Required evidence
# before claiming a novel isoform.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
ISOQ_DIR="${BASE_DIR}/data/isoquant_targeted"
OUT_DIR="${BASE_DIR}/data/sqanti3"
mkdir -p "${OUT_DIR}" "$(pwd)/logs"

MM10_REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/mm10.fa"
GENCODE_GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/gencode.vM25.annotation.gtf"

# IsoQuant writes the novel transcriptome GTF here:
NOVEL_GTF="${ISOQ_DIR}/OUT.transcript_models.gtf"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

sqanti3_qc.py \
    "${NOVEL_GTF}" \
    "${GENCODE_GTF}" \
    "${MM10_REF}" \
    --output ribostamp_srrm3 \
    --dir "${OUT_DIR}" \
    --cpus "${SLURM_NTASKS}" \
    --report skip

echo "=== Done: $(date) ==="
echo "Key outputs:"
echo "  ${OUT_DIR}/ribostamp_srrm3_classification.txt   # per-isoform class + artifact flags"
echo "  ${OUT_DIR}/ribostamp_srrm3_corrected.gtf        # cleaned novel transcriptome"
