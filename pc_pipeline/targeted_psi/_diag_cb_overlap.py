#!/usr/bin/env python3
"""Diagnostic: do long-read and short-read libraries share cell barcodes?

If they do (both run on the same 10x emulsion → same droplets), we'd expect
substantial overlap between the two whitelists. If they don't, the libraries
were run on separate emulsions and barcodes are not comparable.
"""
import pandas as pd
from pathlib import Path

META_DIR = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata")

MICE = ["mouse1", "mouse2", "mouse3"]
LR_GSMS = {"mouse1": "GSM9380799", "mouse2": "GSM9380800", "mouse3": "GSM9380801"}
SR_GSMS = {"mouse1": "GSM9380796", "mouse2": "GSM9380797", "mouse3": "GSM9380798"}

print("Loading whitelists per mouse...")
print(f"\n{'Mouse':<10} {'LR cells':>10} {'SR cells':>10} {'Intersect':>12} {'LR uniq':>10} {'SR uniq':>10}")
print("-" * 70)

for mouse in MICE:
    lr_path = list(META_DIR.glob(f"{LR_GSMS[mouse]}_longread_normed_counts*_obs.tsv"))[0]
    sr_path = list(META_DIR.glob(f"{SR_GSMS[mouse]}_shortread_*_obs.tsv"))[0]

    lr_obs = pd.read_csv(lr_path, sep="\t")
    sr_obs = pd.read_csv(sr_path, sep="\t")

    # Strip ANY suffix: long-read uses "_sampleN", short-read uses "-1"
    import re as _re
    def strip_suffix(s):
        # Remove first non-ACGT char and everything after — barcodes are pure ACGT
        m = _re.match(r"^([ACGT]+)", str(s))
        return m.group(1) if m else str(s)
    lr_bc = set(lr_obs["barcode"].astype(str).map(strip_suffix))
    sr_bc = set(sr_obs["barcode"].astype(str).map(strip_suffix))

    intersect = lr_bc & sr_bc
    lr_only = lr_bc - sr_bc
    sr_only = sr_bc - lr_bc
    print(f"{mouse:<10} {len(lr_bc):>10,} {len(sr_bc):>10,} {len(intersect):>12,} {len(lr_only):>10,} {len(sr_only):>10,}")

print("\n--- Sample barcodes from each ---")
for mouse in MICE[:1]:  # just show mouse1 examples
    lr_path = list(META_DIR.glob(f"{LR_GSMS[mouse]}_longread_normed_counts*_obs.tsv"))[0]
    sr_path = list(META_DIR.glob(f"{SR_GSMS[mouse]}_shortread_*_obs.tsv"))[0]
    lr_obs = pd.read_csv(lr_path, sep="\t")
    sr_obs = pd.read_csv(sr_path, sep="\t")

    print(f"\n{mouse} long-read barcodes (5 examples):")
    print("  ", lr_obs["barcode"].head(5).tolist())
    print(f"\n{mouse} short-read barcodes (5 examples):")
    print("  ", sr_obs["barcode"].head(5).tolist())

# Check: do the cb_extracted (raw) barcodes from long-read per_read.tsv
# match anything in the short-read whitelist?
print("\n--- Cross-check: long-read cb_extracted vs short-read whitelist ---")
PER_READ = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/per_read/all_samples.per_read.tsv")
SRR_TO_MOUSE = {"SRR36480452": "mouse3", "SRR36480453": "mouse2", "SRR36480454": "mouse1"}

# Build full short-read whitelist per mouse (all cell types, not just progenitors)
sr_whitelist_per_mouse = {}
for mouse in MICE:
    sr_path = list(META_DIR.glob(f"{SR_GSMS[mouse]}_shortread_*_obs.tsv"))[0]
    sr_obs = pd.read_csv(sr_path, sep="\t")
    sr_bcs = set(sr_obs["barcode"].astype(str).str.split("_").str[0])
    sr_whitelist_per_mouse[mouse] = sr_bcs

n_total = n_with_cb = n_match_lr = n_match_sr_only = n_no_match = 0
print("Streaming per_read.tsv (this may take a minute)...")
for chunk in pd.read_csv(PER_READ, sep="\t", chunksize=100_000):
    for _, row in chunk.iterrows():
        n_total += 1
        srr = row["sample"]
        mouse = SRR_TO_MOUSE.get(srr)
        if mouse is None:
            continue
        cb = row["cb_extracted"]
        cb_matched_lr = row["cb_matched"]
        if pd.isna(cb) or cb == "":
            continue
        n_with_cb += 1
        if cb_matched_lr and not pd.isna(cb_matched_lr) and cb_matched_lr != "":
            n_match_lr += 1
        if cb in sr_whitelist_per_mouse[mouse]:
            if not (cb_matched_lr and not pd.isna(cb_matched_lr) and cb_matched_lr != ""):
                n_match_sr_only += 1

print(f"\n  Total per_read rows:                              {n_total:,}")
print(f"  With non-empty cb_extracted:                      {n_with_cb:,}")
print(f"  Matched against LONG-read whitelist (cb_matched): {n_match_lr:,}")
print(f"  Matched against SHORT-read whitelist (NOT LR):    {n_match_sr_only:,}")
