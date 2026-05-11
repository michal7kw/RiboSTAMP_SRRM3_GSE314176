#!/bin/bash
#SBATCH --job-name=ribo_05_isoq_tgt
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/05_isoq_tgt_%j.err
#SBATCH --output=./logs/05_isoq_tgt_%j.out

# =============================================================================
# PRIMARY novel-isoform discovery at the Srrm3 locus.
#
# Strategy: subset all 3 HiFi BAMs to Srrm3 ± 10 kb (hpc/srrm3_locus.bed),
# then run IsoQuant with novel-transcript discovery enabled.
# ~1 h on 16 CPU / 32 GB vs ~18 h full-transcriptome.
#
# Full-transcriptome run (06_isoquant_fulltx.sh) is on-demand only.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
ALIGN_DIR="${BASE_DIR}/data/aligned_mm10"
OUT_DIR="${BASE_DIR}/data/isoquant_targeted"
SUBSET_DIR="${BASE_DIR}/data/srrm3_subset"
mkdir -p "${OUT_DIR}" "${SUBSET_DIR}" "$(pwd)/logs"

MM10_REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/mm10.fa"
# TODO: provision Gencode vM25 GTF (mouse, mm10-matched) at this path:
GENCODE_GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/mm10/gencode.vM25.annotation.gtf"
LOCUS_BED="$(dirname "$0")/srrm3_locus.bed"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

echo "=== Step 1: subset all mm10 BAMs to Srrm3 locus ==="
SUBSET_BAMS=()
for BAM in "${ALIGN_DIR}"/*.mm10.bam; do
  [ -f "${BAM}" ] || continue
  NAME="$(basename "${BAM}" .mm10.bam)"
  OUT_BAM="${SUBSET_DIR}/${NAME}.srrm3.bam"
  if [ ! -s "${OUT_BAM}" ]; then
    samtools view -b -L "${LOCUS_BED}" "${BAM}" > "${OUT_BAM}"
    samtools index "${OUT_BAM}"
  fi
  SUBSET_BAMS+=("${OUT_BAM}")
done
echo "Subsetted ${#SUBSET_BAMS[@]} BAMs."

if [ "${#SUBSET_BAMS[@]}" -eq 0 ]; then
  echo "ERROR: no aligned BAMs in ${ALIGN_DIR}. Run 04_align_pbmm2.sh first." >&2
  exit 1
fi

echo "=== Step 2: IsoQuant with novel discovery enabled (high sensitivity) ==="
# --sqanti_output                       -> in-line transcript classification.
# --read_group tag:CB                   -> per-cell quantification keyed on 10x barcode.
# --model_construction_strategy sensitive_pacbio
#   IsoQuant defaults are conservative — they prefer to assign reads to existing
#   annotated transcripts and only emit a novel transcript with strong evidence.
#   We are explicitly hunting a novel isoform (lab short-read anchor at
#   chr5:135,898,574-135,898,652 mm39, 79 bp), so we want sensitivity. False
#   positives at this stage get filtered by SQANTI3 in step 07.
# --genedb_output
#   Emit an updated GTF that includes both annotated AND novel transcripts.
#   This is what gets fed to SQANTI3 + liftOver downstream.
# --complete_genedb is intentionally OMITTED -> novel transcripts survive.
#   See docs/LEARNING_06_novelty_layered_defense.md.

isoquant \
    --reference "${MM10_REF}" \
    --genedb "${GENCODE_GTF}" \
    --bam "${SUBSET_BAMS[@]}" \
    --data_type pacbio_ccs \
    --model_construction_strategy sensitive_pacbio \
    --genedb_output \
    --output "${OUT_DIR}" \
    --threads "${SLURM_NTASKS}" \
    --sqanti_output \
    --count_exons \
    --read_group tag:CB

echo "=== Done: $(date) ==="
echo "Outputs in ${OUT_DIR}:"
ls -lh "${OUT_DIR}"
