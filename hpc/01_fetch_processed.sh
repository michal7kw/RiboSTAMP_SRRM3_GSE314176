#!/bin/bash
#SBATCH --job-name=ribo_01_processed
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error=./logs/01_processed_%j.err
#SBATCH --output=./logs/01_processed_%j.out

# =============================================================================
# Download GSE314176 processed data.
#
# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314176
# Size:   one 15.4 GB TAR archive — GSE314176_RAW.tar
# Contents (per GEO supplementary table): H5AD, MTX, TSV, TXT.
#
# Sufficient for Deliverable 2 (interneuron isoform list) on its own.
# Also produces obs.tsv needed by Deliverable 1's cell-type mapping.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/RiboSTAMP_SRRM3_GSE314176"
PROC_DIR="${BASE_DIR}/data/processed"
META_DIR="${BASE_DIR}/data/metadata"
mkdir -p "${PROC_DIR}" "${META_DIR}" "$(pwd)/logs"

HPC_SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/longread_iso

TAR_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE314176&format=file"
TAR_FILE="${PROC_DIR}/GSE314176_RAW.tar"

echo "=== Step 1: download 15.4 GB TAR from GEO ==="
if [ ! -f "${TAR_FILE}" ]; then
  # -c = resume, -t 5 = 5 retries, -O = explicit output name (GEO's dynamic
  # endpoint otherwise saves as "file" without an extension).
  wget -c -t 5 --tries=5 -O "${TAR_FILE}" "${TAR_URL}"
else
  echo "TAR already downloaded at ${TAR_FILE}; skipping."
fi

echo "=== Step 2: extract ==="
EXTRACT_DIR="${PROC_DIR}/GSE314176_RAW"
mkdir -p "${EXTRACT_DIR}"
# Only extract if there's nothing there yet — tar -k would also work but is
# less predictable across tar implementations.
if [ -z "$(ls -A "${EXTRACT_DIR}")" ]; then
  tar -xvf "${TAR_FILE}" -C "${EXTRACT_DIR}"
else
  echo "Already extracted to ${EXTRACT_DIR}; skipping."
fi

echo "=== Step 3: inspect what we got ==="
echo "Top-level files:"
ls -lh "${EXTRACT_DIR}"

H5AD_FILES=("${EXTRACT_DIR}"/*.h5ad)
MTX_FILES=("${EXTRACT_DIR}"/*.mtx.gz "${EXTRACT_DIR}"/*.mtx)
CSV_FILES=("${EXTRACT_DIR}"/*.csv.gz "${EXTRACT_DIR}"/*.csv)

echo ""
echo "H5AD files found: ${#H5AD_FILES[@]}"
echo "MTX files found:  ${#MTX_FILES[@]}"
echo "CSV files found:  ${#CSV_FILES[@]}"

echo "=== Step 4: gunzip anything compressed in place ==="
# Leave originals alone; produce plain .mtx / .csv / .tsv the local R scripts expect.
for f in "${EXTRACT_DIR}"/*.gz; do
  [ -f "${f}" ] || continue
  out="${f%.gz}"
  if [ ! -f "${out}" ]; then
    gunzip -k "${f}"
  fi
done

echo "=== Step 5: explode H5AD (if present) ==="
# Prefer H5AD if available — it bundles obs/var/layers in one file.
# explode_h5ad.py normalizes the orientation to features x cells and writes
# plain files the R scripts can readMM() directly.
for h5ad in "${H5AD_FILES[@]}"; do
  [ -f "${h5ad}" ] || continue
  echo "Exploding ${h5ad} ..."
  python "${HPC_SCRIPT_DIR}/explode_h5ad.py" "${h5ad}" "${META_DIR}"
done

if [ ! -f "${META_DIR}/obs.tsv" ]; then
  echo ""
  echo "WARNING: no obs.tsv produced. Either the TAR contained no .h5ad, or"
  echo "the h5ad didn't carry an obs DataFrame. Manually inspect:"
  echo "  ${EXTRACT_DIR}"
  echo "and adapt the local/01_parse_geo_metadata.R input path."
fi

echo "=== Done: $(date) ==="
echo "Processed outputs in ${EXTRACT_DIR}"
echo "Cell metadata + exploded matrices in ${META_DIR}"
