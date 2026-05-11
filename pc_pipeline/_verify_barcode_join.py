#!/usr/bin/env python3
"""Verify shortread/longread barcode overlap after normalizing both to
the bare 16-bp 10x barcode (strip any -N or _sampleN suffix)."""
import re
import pandas as pd

M = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata"

PAIRS = [("GSM9380796", "GSM9380799", "mouse1"),
         ("GSM9380797", "GSM9380800", "mouse2"),
         ("GSM9380798", "GSM9380801", "mouse3")]

# Strip "-1", "_sample1" etc. — keep just the bare 16-mer
SUFFIX_RE = re.compile(r"[_\-].*$")
def strip_suffix(b: str) -> str:
    return SUFFIX_RE.sub("", b)


print(f"{'mouse':6s}  {'sr':>5s}  {'lr':>5s}  {'shared':>6s}  {'%lr':>5s}  rescued_OPC  rescued_iDG")
print("-" * 70)
for sr_gsm, lr_gsm, mname in PAIRS:
    sr = pd.read_csv(f"{M}/{sr_gsm}_shortread_epr_adata_{mname}__obs.tsv", sep="\t")
    lr = pd.read_csv(f"{M}/{lr_gsm}_longread_normed_counts_transcript_adata_{mname}__obs.tsv", sep="\t")

    sr["bare"] = sr["barcode"].astype(str).map(strip_suffix)
    lr["bare"] = lr["barcode"].astype(str).map(strip_suffix)

    sr_bc = set(sr["bare"])
    lr_bc = set(lr["bare"])
    common = sr_bc & lr_bc

    sr_in_common = sr[sr["bare"].isin(common)]
    cag = "Cell Assignments Grouped"
    n_opc = (sr_in_common[cag] == "Glia - OPC").sum()
    n_idg = (sr_in_common[cag] == "Neuron - DG - Immature").sum()

    pct = 100 * len(common) / max(len(lr_bc), 1)
    print(f"{mname:6s}  {len(sr_bc):>5d}  {len(lr_bc):>5d}  {len(common):>6d}  {pct:>4.1f}%  {n_opc:>11d}  {n_idg:>11d}")

# Sanity: a few examples of bare barcodes that match across both
print()
print("=== mouse1 — sample shared bare barcodes ===")
sr = pd.read_csv(f"{M}/GSM9380796_shortread_epr_adata_mouse1__obs.tsv", sep="\t")
lr = pd.read_csv(f"{M}/GSM9380799_longread_normed_counts_transcript_adata_mouse1__obs.tsv", sep="\t")
sr["bare"] = sr["barcode"].map(strip_suffix)
lr["bare"] = lr["barcode"].map(strip_suffix)
sr["src"] = "sr"
lr["src"] = "lr"
common = set(sr["bare"]) & set(lr["bare"])
sample = list(common)[:5]
for bc in sample:
    sr_row = sr[sr["bare"] == bc].iloc[0]
    lr_row = lr[lr["bare"] == bc].iloc[0]
    print(f"  {bc}:  sr={sr_row['barcode']!s} ({sr_row['Cell Assignments Grouped']!s}) | lr={lr_row['barcode']!s} ({lr_row['Cell Assignments Grouped']!s})")
