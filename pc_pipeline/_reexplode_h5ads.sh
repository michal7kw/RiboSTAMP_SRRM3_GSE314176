#!/bin/bash
# One-off helper — wipe stale single-h5ad outputs in metadata/ then re-explode
# the six long-read h5ads (3 mice × {normed_counts, EditsC}).
set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

META=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata
HPC=/mnt/d/Github/SRF/Elisa/RiboSTAMP_SRRM3_GSE314176/hpc
EXTRACT=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/processed/GSE314176_RAW

echo "=== clean stale unnamespaced outputs ==="
rm -fv "$META"/obs.tsv "$META"/var.tsv "$META"/counts.mtx "$META"/layer_counts.mtx 2>/dev/null || true

echo ""
echo "=== re-explode the 6 long-read h5ads (Linda Track A + EditsC) ==="
python "$HPC/explode_h5ad.py" "$META" \
    "$EXTRACT/GSM9380799_longread_normed_counts_transcript_adata_mouse1.h5ad" \
    "$EXTRACT/GSM9380800_longread_normed_counts_transcript_adata_mouse2.h5ad" \
    "$EXTRACT/GSM9380801_longread_normed_counts_transcript_adata_mouse3.h5ad" \
    "$EXTRACT/GSM9380799_longread_EditsC_transcript_adata_mouse1.h5ad" \
    "$EXTRACT/GSM9380800_longread_EditsC_transcript_adata_mouse2.h5ad" \
    "$EXTRACT/GSM9380801_longread_EditsC_transcript_adata_mouse3.h5ad"

echo ""
echo "=== explode the 3 SHORT-READ h5ads (for OPC+iDG cell-type rescue, Track B) ==="
# These contain richer cell-type annotations including OPC and iDG-Immature
# that the long-read h5ads filtered out. Same biological mice, paired by
# mouse1/2/3 -> we can join on barcode to recover progenitor labels.
python "$HPC/explode_h5ad.py" "$META" \
    "$EXTRACT/GSM9380796_shortread_epr_adata_mouse1.h5ad" \
    "$EXTRACT/GSM9380797_shortread_epr_adata_mouse2.h5ad" \
    "$EXTRACT/GSM9380798_shortread_epr_adata_mouse3.h5ad"

echo ""
echo "=== metadata dir contents ==="
ls -lh "$META"
