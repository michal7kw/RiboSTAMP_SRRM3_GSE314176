#!/bin/bash
# Quick check on parallel step-03 progress.
source /home/michal/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate longread_iso 2>/dev/null

BASE=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi
for srr in SRR36480452 SRR36480453 SRR36480454; do
    total=$(samtools view -c "${BASE}/aligned_locus/${srr}.locus.bam" 2>/dev/null)
    processed=$(wc -l < "${BASE}/per_read/${srr}.per_read.tsv" 2>/dev/null)
    echo "${srr}: total=${total}  processed=${processed}"
done

echo
echo "--- merged TSV (if any) ---"
ls -lh "${BASE}/per_read/all_samples.per_read.tsv" 2>/dev/null || echo "not yet merged"

echo
echo "--- results dir ---"
ls -lh "${BASE}/results/" 2>/dev/null || echo "results not yet present"
