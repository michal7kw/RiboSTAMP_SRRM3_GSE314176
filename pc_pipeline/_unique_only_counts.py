#!/usr/bin/env python3
"""
Tighter version of the per-cluster expression analysis using ONLY uniquely
assigned reads (assignment_type in {"unique", "unique_minor_difference"}).

This addresses the previous-turn caveat that re-quant counts overlap between
similar isoforms — the unique-only counts cannot be double-assigned, so each
read contributes to exactly one transcript.
"""
from __future__ import annotations

import gzip
from pathlib import Path
import pandas as pd

BASE = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data")
ISOQ_OUT = BASE / "isoquant_targeted_quant" / "OUT"
TGT_PERREAD = BASE / "targeted_psi" / "per_read" / "all_samples.per_read.tsv"
META = BASE / "metadata"
RES = BASE / "isoquant_targeted" / "results"

SRR_TO_GSM = {
    "SRR36480452": "GSM9380801",
    "SRR36480453": "GSM9380800",
    "SRR36480454": "GSM9380799",
}

TRANSCRIPTS = {
    "ENSMUST00000005077.6":    "Srrm3-201 (canonical)",
    "ENSMUST00000144211.1":    "Srrm3-204 (ref of t185)",
    "transcript185.chr5.nnic": "novel t185 (skips Ex15, non-coding)",
    "transcript11.chr5.nnic":  "novel t11 (11 exons, NMD)",
    "transcript157.chr5.nnic": "novel t157 (extra novel intron)",
}

# Load read_assignments (only header + isoform/assignment_type)
print("Loading read_assignments...")
ra = []
with gzip.open(ISOQ_OUT / "OUT.read_assignments.tsv.gz", "rt") as fh:
    header = None
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            ri = header.index("read_id")
            ii = header.index("isoform_id")
            ai = header.index("assignment_type")
            continue
        ra.append((parts[ri], parts[ii], parts[ai]))
ra = pd.DataFrame(ra, columns=["read_id", "isoform_id", "assignment_type"])

# Filter to unique assignments only
unique_only = ra[ra["assignment_type"].isin(["unique", "unique_minor_difference"])].copy()
print(f"  total assignments: {len(ra):,}")
print(f"  unique-only assignments: {len(unique_only):,}")

# Per-transcript read totals (unique-only)
print("\n=== Unique-only read counts per transcript ===")
counts_per_iso = unique_only["isoform_id"].value_counts()
for tid, label in TRANSCRIPTS.items():
    n = counts_per_iso.get(tid, 0)
    print(f"  {tid:30s}  {n:>6d}  ({label})")

# Subset and join
unique_subset = unique_only[unique_only["isoform_id"].isin(TRANSCRIPTS.keys())].copy()
print(f"\n  unique-only reads on transcripts of interest: {len(unique_subset):,}")

# Join CB
pr = pd.read_csv(TGT_PERREAD, sep="\t",
                 dtype={"cb_matched": "string", "sample": "string"})
pr["cb_matched"] = pr["cb_matched"].fillna("")
pr = pr[["sample", "read_id", "cb_matched"]].rename(columns={"sample": "srr"})

joined = unique_subset.merge(pr, on="read_id", how="left")
cb_recov = (joined["cb_matched"].fillna("") != "").sum()
print(f"  unique-only reads with CB: {cb_recov:,}/{len(joined):,} ({100*cb_recov/max(1,len(joined)):.1f}%)")

# Join cell type
obs_parts = []
for srr, gsm in SRR_TO_GSM.items():
    p = next(META.glob(f"{gsm}_longread_normed_counts_*__obs.tsv"), None)
    df = pd.read_csv(p, sep="\t")
    df["srr"] = srr
    df["bare_barcode"] = df["barcode"].astype(str).str.split("_", n=1).str[0]
    df = df.rename(columns={"Cell Assignments Grouped": "cell_type"})
    obs_parts.append(df[["srr", "bare_barcode", "cell_type"]])
obs = pd.concat(obs_parts, ignore_index=True)

joined = joined.merge(
    obs, left_on=["srr", "cb_matched"],
    right_on=["srr", "bare_barcode"], how="left"
)

# Filter to (CB-matched + cell-type-assigned) reads
kept = joined[(joined["cb_matched"].fillna("") != "") &
              joined["cell_type"].notna()]
print(f"  unique reads with both CB and cell-type: {len(kept):,}")

# Pivot table
pivot = (kept.groupby(["cell_type", "isoform_id"]).size()
              .unstack(fill_value=0))
total_per_cluster = pivot.sum(axis=1)
pivot_pct = pivot.div(total_per_cluster, axis=0) * 100

# Order columns
col_order = [c for c in [
    "ENSMUST00000005077.6", "ENSMUST00000144211.1",
    "transcript185.chr5.nnic", "transcript11.chr5.nnic", "transcript157.chr5.nnic",
] if c in pivot.columns]
pivot = pivot[col_order]
pivot_pct = pivot_pct[col_order]
pivot["TOTAL_unique"] = total_per_cluster

# Order rows
row_order = ["BBB - Endo", "Glia - Astro", "Glia - Oligo",
             "Neuron - DG", "Neuron - CA1", "Neuron - CA3", "Neuron - GABA"]
row_order = [r for r in row_order if r in pivot.index]
pivot = pivot.reindex(row_order)
pivot_pct = pivot_pct.reindex(row_order)

print("\n=== UNIQUE-only counts per cluster (no double-assignments) ===")
print(pivot.to_string())
print("\n=== UNIQUE-only % composition ===")
print(pivot_pct.round(2).to_string())

pivot.to_csv(RES / "novel_isoforms_per_cluster_UNIQUE_counts.tsv", sep="\t")
pivot_pct.to_csv(RES / "novel_isoforms_per_cluster_UNIQUE_pct.tsv", sep="\t")
print("\nWrote unique-only outputs.")
