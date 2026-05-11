#!/usr/bin/env python3
"""For one author barcode, search its forward + RC sequence in actual reads
to verify what orientation/position it actually appears at."""
import re
import pandas as pd
import pysam

BAM = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/aligned_locus/SRR36480452.locus.bam"
WHITELIST = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata/GSM9380801_longread_normed_counts_transcript_adata_mouse3__obs.tsv"

REVCOMP = str.maketrans("ACGTNacgtn", "TGCANTGCAN")
def revcomp(s): return s.translate(REVCOMP)[::-1]

# Load all author barcodes
wl_df = pd.read_csv(WHITELIST, sep="\t")
barcodes = sorted(set(wl_df["barcode"].astype(str).str.split("_", n=1).str[0]))
print(f"Total barcodes: {len(barcodes)}")
print(f"Sample of 5: {barcodes[:5]}")

# Build a regex of all barcodes (forward + RC) for fast multi-match
fwd_set = set(barcodes)
rc_set = set(revcomp(b) for b in barcodes)
print(f"Forward set: {len(fwd_set)}, RC set: {len(rc_set)}")

# Test a SINGLE barcode against many reads
target_fwd = barcodes[0]   # first author barcode
target_rc = revcomp(target_fwd)
print(f"\nSearching for {target_fwd} (forward) and {target_rc} (RC) in BAM reads...")

n_total = 0
n_fwd_hit = 0
n_rc_hit = 0
fwd_examples = []
rc_examples = []

with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        n_total += 1
        seq = read.query_sequence
        if seq is None: continue
        # Try as-is
        if target_fwd in seq:
            n_fwd_hit += 1
            if len(fwd_examples) < 3:
                idx = seq.index(target_fwd)
                ctx = seq[max(0, idx-30):idx + 16 + 30]
                fwd_examples.append(ctx)
        if target_rc in seq:
            n_rc_hit += 1
            if len(rc_examples) < 3:
                idx = seq.index(target_rc)
                ctx = seq[max(0, idx-30):idx + 16 + 30]
                rc_examples.append(ctx)
        if n_total >= 50000: break

print(f"\nFrom {n_total} reads:")
print(f"  exact match {target_fwd!r} (fwd): {n_fwd_hit}")
print(f"  exact match {target_rc!r} (rc): {n_rc_hit}")

print(f"\nForward hits with context (30bp flanks):")
for ex in fwd_examples:
    print(f"  ...{ex}...")
print(f"\nRC hits with context (30bp flanks):")
for ex in rc_examples:
    print(f"  ...{ex}...")

# How often does ANY author barcode (fwd or rc) appear in reads?
print(f"\n=== batch search: any of {len(barcodes)} barcodes in first 5000 reads ===")
import time
start = time.time()
n_total = 0
n_any_match = 0
n_fwd_match = 0
n_rc_match = 0
with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        n_total += 1
        seq = read.query_sequence
        if seq is None: continue
        # 16-mer scan: walk every 16-bp window, check each window against fwd or rc set
        # That's expensive; instead extract candidates near polyA/T like the script does
        # ...for now just search a sample of barcodes
        for i in range(0, len(seq) - 16):
            window = seq[i:i + 16]
            if window in fwd_set or window in rc_set:
                n_any_match += 1
                if window in fwd_set: n_fwd_match += 1
                if window in rc_set: n_rc_match += 1
                break
        if n_total >= 5000: break

dur = time.time() - start
print(f"  in {n_total} reads ({dur:.1f}s): {n_any_match} have exact match to any author barcode")
print(f"    fwd matches: {n_fwd_match}")
print(f"    rc matches:  {n_rc_match}")
