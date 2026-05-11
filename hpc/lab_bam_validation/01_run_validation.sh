#!/bin/bash
#SBATCH --job-name=ribo_lab_bam_psi
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/lab_bam_psi_%j.err
#SBATCH --output=./logs/lab_bam_psi_%j.out

# =============================================================================
# Validate the junction-PSI calculator against the lab's mm39 N2A bulk BAMs.
#
# Pass criterion: per-group PSI matches the published anchor
#   Parental ~57%   Neg ~40%   Pos ~27%   KO ~5%
# (from 90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md)
#
# If recovered, the calculator is validated for cross-context use (Surface A,
# Surface B). See docs/HPC_NEXT_STEPS.md Task 1 for the wider context.
#
# Note on log output: 03_count_junctions.py prints "Junctions (mm10):" — that
# header is hardcoded but the actual coordinates come from the SR_* env vars
# below, which are mm39. Cosmetic only.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/RiboSTAMP_SRRM3_GSE314176"
LAB_BAMS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/results/02_aligned"
OUT_DIR="${BASE_DIR}/data/lab_bam_validation"
COUNTER="${BASE_DIR}/pc_pipeline/sr_junction_psi/03_count_junctions.py"
mkdir -p "${OUT_DIR}" "$(pwd)/logs"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

# mm39 cassette coordinates — verified against the lab's anchor analysis report
export SR_CHROM="chr5"
export SR_INC1_LEFT=135898148
export SR_INC1_RIGHT=135898574
export SR_INC2_LEFT=135898652
export SR_INC2_RIGHT=135901934
export SR_JUNC_TOL=5

declare -A SAMPLE_TO_GROUP=(
    [1]="Parental"  [2]="Parental"  [3]="Parental"
    [4]="Neg"       [5]="Neg"       [6]="Neg"
    [7]="Pos"       [8]="Pos"       [9]="Pos"
    [13]="KO"       [14]="KO"       [15]="KO"
)

PER_SAMPLE_TSV="${OUT_DIR}/lab_n2a_psi.tsv"
: > "${PER_SAMPLE_TSV}"
HEADER_WRITTEN=0

echo "============================================"
echo "Lab N2A bulk PSI validation"
echo "Start:       $(date)"
echo "BAMs dir:    ${LAB_BAMS_DIR}"
echo "Out dir:     ${OUT_DIR}"
echo "Calculator:  ${COUNTER}"
echo "Coords mm39: ${SR_CHROM}:${SR_INC1_LEFT}-${SR_INC2_RIGHT} (cassette ${SR_INC1_RIGHT}-${SR_INC2_LEFT})"
echo "============================================"

for SAMPLE in 1 2 3 4 5 6 7 8 9 13 14 15; do
    GROUP="${SAMPLE_TO_GROUP[$SAMPLE]}"
    BAM="${LAB_BAMS_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    OUT="${OUT_DIR}/sample_${SAMPLE}_${GROUP}_psi.tsv"

    echo ""
    echo "=== Sample ${SAMPLE} (${GROUP}) ==="

    if [[ ! -s "${BAM}" ]]; then
        echo "  [skip] BAM not found or empty: ${BAM}"
        continue
    fi

    python "${COUNTER}" \
        --bam "${BAM}" \
        --output "${OUT}" \
        --sample-name "sample_${SAMPLE}_${GROUP}"

    if [[ ${HEADER_WRITTEN} -eq 0 ]]; then
        cat "${OUT}" >> "${PER_SAMPLE_TSV}"
        HEADER_WRITTEN=1
    else
        tail -n +2 "${OUT}" >> "${PER_SAMPLE_TSV}"
    fi
done

echo ""
echo "============================================"
echo "Per-group PSI aggregation"
echo "============================================"

export PER_SAMPLE_TSV
export OUT_DIR

python <<'PYEOF'
import os
import pandas as pd

per_sample_tsv = os.environ["PER_SAMPLE_TSV"]
out_dir = os.environ["OUT_DIR"]

df = pd.read_csv(per_sample_tsv, sep="\t")
df["group"] = df["sample"].str.extract(r"_(Parental|Neg|Pos|KO)$")[0]

agg = (
    df.groupby("group", sort=False)
      .agg(
          n_samples=("sample", "count"),
          n_inc1=("n_inc1", "sum"),
          n_inc2=("n_inc2", "sum"),
          n_skip=("n_skip", "sum"),
      )
      .reset_index()
)
agg["psi"] = (agg["n_inc1"] + agg["n_inc2"]) / (
    agg["n_inc1"] + agg["n_inc2"] + 2 * agg["n_skip"]
)
agg["psi_pct"] = (agg["psi"] * 100).round(1)

# Reorder rows to match the lab's reporting order
order = ["Parental", "Neg", "Pos", "KO"]
agg["__ord"] = agg["group"].map({g: i for i, g in enumerate(order)})
agg = agg.sort_values("__ord").drop(columns="__ord").reset_index(drop=True)

print()
print(agg.to_string(index=False))
print()
print("Lab anchor (90-1239779069):  Parental=57%  Neg=40%  Pos=27%  KO=5%")
print("Pass bands:                  Parental 30-75%  Neg 25-50%  Pos 18-37%  KO 0-15%")
print()

out_path = os.path.join(out_dir, "lab_n2a_psi_per_group.tsv")
agg.to_csv(out_path, sep="\t", index=False)
print(f"Per-group table written to: {out_path}")
PYEOF

echo ""
echo "=== Done: $(date) ==="
