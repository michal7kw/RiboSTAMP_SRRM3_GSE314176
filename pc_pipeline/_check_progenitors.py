#!/usr/bin/env python3
"""Search both long-read AND short-read obs.tsv files for any cell labels
that look like progenitors / immature / proliferating populations.

Specifically searching for: OPC, NPC, iDG, immature, progenitor, proliferating,
stem, dividing, neuroblast, pre-, _imm, COP (committed oligo precursor).
"""
from pathlib import Path
from collections import Counter
import re

import pandas as pd

META_DIR = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata")

LR_FILES = {
    "mouse1 LR": META_DIR / "GSM9380799_longread_normed_counts_transcript_adata_mouse1__obs.tsv",
    "mouse2 LR": META_DIR / "GSM9380800_longread_normed_counts_transcript_adata_mouse2__obs.tsv",
    "mouse3 LR": META_DIR / "GSM9380801_longread_normed_counts_transcript_adata_mouse3__obs.tsv",
}
SR_FILES = {
    "mouse1 SR": META_DIR / "GSM9380796_shortread_epr_adata_mouse1__obs.tsv",
    "mouse2 SR": META_DIR / "GSM9380797_shortread_epr_adata_mouse2__obs.tsv",
    "mouse3 SR": META_DIR / "GSM9380798_shortread_epr_adata_mouse3__obs.tsv",
}

PROG_HINTS = re.compile(
    r"\b(opc|cop|npc|idg|immature|imm$|imm_|_imm|progenitor|prolif|stem|"
    r"dividing|neuroblast|pre-|precursor|cycling|nsc)\b",
    re.IGNORECASE,
)

def scan_file(label, path):
    if not path.exists():
        print(f"   [missing] {label}: {path}")
        return None, None
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"\n--- {label}  (n={len(df):,} cells; columns: {list(df.columns)}) ---")
    grouped = df.get("Cell Assignments Grouped")
    sub = df.get("Cell Assignment")
    if grouped is not None:
        print(f"  unique 'Cell Assignments Grouped': {sorted(grouped.dropna().unique())}")
    if sub is not None:
        print(f"  unique 'Cell Assignment' (subtypes, {sub.nunique()} total):")
        for x in sorted(sub.dropna().unique()):
            tag = " <<< PROGENITOR-LIKE" if PROG_HINTS.search(str(x)) else ""
            n = (sub == x).sum()
            print(f"      {str(x):<35s}{n:>6d}{tag}")
    return grouped, sub

print("=" * 70)
print("LONG-READ obs files (used for our PSI)")
print("=" * 70)
for label, path in LR_FILES.items():
    scan_file(label, path)

print("\n" + "=" * 70)
print("SHORT-READ obs files (the *parallel* sn/scRNA-seq from the same paper)")
print("=" * 70)
sr_grouped_pool = Counter()
sr_sub_pool = Counter()
for label, path in SR_FILES.items():
    grouped, sub = scan_file(label, path)
    if grouped is not None:
        sr_grouped_pool.update(grouped.dropna())
    if sub is not None:
        sr_sub_pool.update(sub.dropna())

print("\n" + "=" * 70)
print("Short-read POOLED — Cell Assignments Grouped")
print("=" * 70)
total = sum(sr_grouped_pool.values())
for ct, n in sorted(sr_grouped_pool.items(), key=lambda x: -x[1]):
    tag = " <<< PROGENITOR-LIKE" if PROG_HINTS.search(str(ct)) else ""
    print(f"   {ct:<35s}{n:>6d}   ({100*n/total:5.1f}%){tag}")
print(f"   {'TOTAL':<35s}{total:>6d}")

print("\n" + "=" * 70)
print("Short-read POOLED — Cell Assignment subtypes (full list)")
print("=" * 70)
total_sub = sum(sr_sub_pool.values())
for ct, n in sorted(sr_sub_pool.items(), key=lambda x: -x[1]):
    tag = " <<< PROGENITOR-LIKE" if PROG_HINTS.search(str(ct)) else ""
    print(f"   {str(ct):<35s}{n:>6d}{tag}")
print(f"   {'TOTAL':<35s}{total_sub:>6d}")
