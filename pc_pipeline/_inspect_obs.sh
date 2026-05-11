#!/bin/bash
set -euo pipefail
META=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata
for f in "$META"/GSM93807{99,9380800,9380801}_longread_normed_counts_transcript_adata_mouse*_obs.tsv; do :; done
# (just list what we have)
ls "$META" | grep "longread_normed.*obs.tsv"

OBS1="$META/GSM9380799_longread_normed_counts_transcript_adata_mouse1__obs.tsv"
echo ""
echo "=== Mouse1 obs.tsv header ==="
head -1 "$OBS1" | tr '\t' '\n' | nl
echo ""
echo "=== first row ==="
sed -n '2p' "$OBS1" | tr '\t' '\n' | nl
echo ""
echo "=== rows ==="
wc -l "$OBS1"
