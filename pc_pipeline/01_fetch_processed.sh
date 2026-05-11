#!/bin/bash
# =============================================================================
# pc_pipeline/01 — download GSE314176_RAW.tar (15.4 GB), extract, explode h5ad.
# Same logic as hpc/01_fetch_processed.sh, no SBATCH, runs in WSL foreground.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

log_step "Step 1 — download GSE314176_RAW.tar (15.4 GB)"
TAR_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE314176&format=file"
TAR_FILE="${PROC_DIR}/GSE314176_RAW.tar"

if [ ! -f "${TAR_FILE}" ]; then
  wget -c -t 5 -O "${TAR_FILE}" "${TAR_URL}"
else
  echo "TAR already present (${TAR_FILE}); skipping download."
fi

log_step "Step 2 — extract"
EXTRACT_DIR="${PROC_DIR}/GSE314176_RAW"
mkdir -p "${EXTRACT_DIR}"
if [ -z "$(ls -A "${EXTRACT_DIR}" 2>/dev/null)" ]; then
  tar -xvf "${TAR_FILE}" -C "${EXTRACT_DIR}"
else
  echo "Already extracted (${EXTRACT_DIR}); skipping."
fi

log_step "Step 3 — gunzip anything compressed"
for f in "${EXTRACT_DIR}"/*.gz; do
  [ -f "${f}" ] || continue
  out="${f%.gz}"
  [ -f "${out}" ] || gunzip -k "${f}"
done

log_step "Step 4 — explode the LONG-READ h5ads we actually need"
HPC_DIR="$(dirname "$(readlink -f "$0")")/../hpc"

# GSE314176 ships many h5ads:
#   *_shortread_epr_adata_mouseN.h5ad             — Illumina/AVITI scRNA, NOT what we use
#   *_longread_normed_counts_transcript_adata_*   — long-read TRANSCRIPTION (Linda Track A)
#   *_longread_EditsC_transcript_adata_*          — long-read TRANSLATION (EditsC bonus column)
#   *_longread_normed_edits_adata_*               — long-read normalized edit RATES (not used now)
#
# We only explode the two we need per mouse. Each h5ad's outputs are
# namespaced by basename (see hpc/explode_h5ad.py) so they don't clobber.

shopt -s nullglob
COUNTS_FILES=("${EXTRACT_DIR}"/*longread_normed_counts_transcript*.h5ad)
EDITSC_FILES=("${EXTRACT_DIR}"/*longread_EditsC_transcript*.h5ad)
shopt -u nullglob

if [ "${#COUNTS_FILES[@]}" -eq 0 ] && [ "${#EDITSC_FILES[@]}" -eq 0 ]; then
  echo "WARNING: no long-read transcript h5ads found under ${EXTRACT_DIR}." >&2
  echo "         Inspect the archive contents and adjust the glob patterns above." >&2
else
  for h5ad in "${COUNTS_FILES[@]}" "${EDITSC_FILES[@]}"; do
    echo ""
    echo "Exploding ${h5ad}"
    python "${HPC_DIR}/explode_h5ad.py" "${META_DIR}" "${h5ad}"
  done
fi

log_step "Done"
ls -lh "${META_DIR}"
