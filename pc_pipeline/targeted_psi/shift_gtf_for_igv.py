import sys
import os

gtf_in = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/reference/mm10/gencode.vM25.annotation.gtf"
gtf_out = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/reference/srrm3_locus_annotation.gtf"
locus_chr = "chr5"
locus_start = 135839720
locus_end = 135899798
seqname_out = "chr5:135839721-135899798"

count = 0
print(f"Reading from {gtf_in}...")
with open(gtf_in, "r") as fin, open(gtf_out, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        
        chrom = parts[0]
        if chrom != locus_chr:
            continue
            
        start = int(parts[3])
        end = int(parts[4])
        
        # Check for overlap
        if start <= locus_end and end > locus_start:
            # Shift coordinates
            new_start = start - locus_start
            new_end = end - locus_start
            
            # Clip to locus boundaries
            if new_start < 1: new_start = 1
            if new_end > (locus_end - locus_start): new_end = (locus_end - locus_start)
            
            parts[0] = seqname_out
            parts[3] = str(new_start)
            parts[4] = str(new_end)
            
            fout.write("\t".join(parts))
            count += 1

print(f"Success! Wrote {count} GTF features to {gtf_out}")
