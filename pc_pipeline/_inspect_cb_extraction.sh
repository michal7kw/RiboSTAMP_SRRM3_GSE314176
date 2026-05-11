#!/bin/bash
# Inspect what CBs we're extracting and why they don't match the whitelist.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

PERREAD=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/per_read/SRR36480452.per_read.tsv
META=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata
WHITELIST_OBS="${META}/GSM9380801_longread_normed_counts_transcript_adata_mouse3__obs.tsv"

echo "=== n cells in author whitelist (mouse3) ==="
awk -F'\t' 'NR>1 {n++} END {print n}' "$WHITELIST_OBS"

echo ""
echo "=== sample 10 author barcodes ==="
awk -F'\t' 'NR>1 && NR<=11 {print $1}' "$WHITELIST_OBS"

echo ""
echo "=== sample 10 extracted candidates from per_read.tsv (where cb_extracted is non-empty) ==="
awk -F'\t' 'NR>1 && $4 != "" {print $1, $4, $5, $6}' "$PERREAD" | head -10

echo ""
echo "=== count rows where cb_extracted is non-empty (CB candidate found) vs total ==="
awk -F'\t' 'NR>1 {n++; if ($4 != "") got++; if ($5 != "") matched++} END {
  printf "  total rows: %d\n", n
  printf "  cb_extracted: %d (%.2f%%)\n", got, 100*got/n
  printf "  cb_matched:   %d (%.2f%%)\n", matched, 100*matched/n
}' "$PERREAD"
