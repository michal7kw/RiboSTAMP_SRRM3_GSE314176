#!/bin/bash
# Launch step 03 (CIGAR + barcode classification) in PARALLEL for 3 SRRs,
# then merge into all_samples.per_read.tsv, then run step 04.
# ~3× faster than the sequential run_all path.
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "${HERE}/00_config.sh"
activate_env

LOG_DIR="${TGT_LOG_DIR}"
mkdir -p "${LOG_DIR}"

ts=$(date +%Y%m%d_%H%M%S)

echo "=== launching 3 parallel step 03 processes ==="
for SRR in SRR36480452 SRR36480453 SRR36480454; do
  log="${LOG_DIR}/03_${SRR}_${ts}.log"
  echo "  ${SRR} -> ${log}"
  python "${HERE}/03_classify_and_barcode.py" --srr "${SRR}" > "${log}" 2>&1 &
done

echo "Waiting for all 3 to finish..."
wait
echo "All step 03 processes done."

echo ""
echo "=== merging individual per_read.tsv into all_samples.per_read.tsv ==="
PERREAD="${TGT_PERREAD_DIR}"
POOLED="${PERREAD}/all_samples.per_read.tsv"

head -1 "${PERREAD}/SRR36480452.per_read.tsv" > "${POOLED}"
for SRR in SRR36480452 SRR36480453 SRR36480454; do
  tail -n +2 "${PERREAD}/${SRR}.per_read.tsv" >> "${POOLED}"
done

n_rows=$(wc -l < "${POOLED}")
echo "Pooled per-read TSV: ${POOLED} ($((n_rows - 1)) rows)"

echo ""
echo "=== running step 04 ==="
python "${HERE}/04_per_cluster_psi.py" 2>&1 | tee "${LOG_DIR}/04_${ts}.log"

echo ""
echo "=== done ==="
ls -lh "${TGT_RESULTS_DIR}/" 2>/dev/null
