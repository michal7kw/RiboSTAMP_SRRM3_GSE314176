#!/usr/bin/env python3
"""
Normalized Srrm3 isoform composition per cell type.

For each cell type, what % of its total Srrm3 reads come from each isoform?
This controls for differing cluster sizes (Glia-Astro has the most cells, so
absolute counts will be high there for *any* isoform).

Reads from the per-read joined table — only reads with a CB and a cell-type
assignment are considered ("UNASSIGNED" reads are excluded).
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

BASE = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/results")
PER_READ = BASE / "novel_isoforms_per_read.tsv"

LABELS = {
    "ENSMUST00000005077.6":    "Srrm3-201 (canonical ref)",
    "ENSMUST00000144211.1":    "Srrm3-204 (ref of t185)",
    "transcript11.chr5.nnic":  "novel t11 (11 exons)",
    "transcript157.chr5.nnic": "novel t157 (extra novel intron)",
    "transcript185.chr5.nnic": "novel t185 (skips Ex15)",
}

df = pd.read_csv(PER_READ, sep="\t", dtype={"cb_matched": "string", "cell_type": "string"})

# only reads with a CB AND author-assigned cell_type
df = df[df["cb_matched"].fillna("") != ""]
df = df[df["cell_type"].notna() & (df["cell_type"] != "UNASSIGNED")]
print(f"Reads with both CB and cell-type assignment: {len(df):,}")

ct_iso = df.groupby(["cell_type", "isoform_id"]).size().unstack(fill_value=0)
ct_total = ct_iso.sum(axis=1)
ct_pct = ct_iso.div(ct_total, axis=0) * 100

# Order columns: canonical first, then novel
col_order = [c for c in [
    "ENSMUST00000005077.6", "ENSMUST00000144211.1",
    "transcript185.chr5.nnic", "transcript11.chr5.nnic", "transcript157.chr5.nnic",
] if c in ct_iso.columns]
ct_iso = ct_iso[col_order]
ct_pct = ct_pct[col_order]
ct_pct["TOTAL_READS"] = ct_total

# Number of cells per cluster (from author obs)
import os
META = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata")
SRR_TO_GSM = {
    "SRR36480452": "GSM9380801",
    "SRR36480453": "GSM9380800",
    "SRR36480454": "GSM9380799",
}
obs_parts = []
for srr, gsm in SRR_TO_GSM.items():
    p = next(META.glob(f"{gsm}_longread_normed_counts_*__obs.tsv"), None)
    obs_parts.append(pd.read_csv(p, sep="\t"))
obs_all = pd.concat(obs_parts, ignore_index=True).rename(columns={
    "Cell Assignments Grouped": "cell_type"
})
n_cells = obs_all.groupby("cell_type").size().rename("n_cells")
ct_pct = ct_pct.join(n_cells, how="left")

print("\n=== Cluster-normalized Srrm3 isoform composition (% of cluster's total Srrm3 reads) ===")
print(ct_pct.round(2).to_string())

print("\n=== Read counts per cluster per isoform ===")
print(ct_iso.to_string())

ct_pct.to_csv(BASE / "novel_isoforms_per_cluster_pct.tsv", sep="\t")
ct_iso.to_csv(BASE / "novel_isoforms_per_cluster_counts.tsv", sep="\t")

# Reads per 100 cells per cluster (a cell-density-normalized read count)
ct_per_100c = ct_iso.div(n_cells, axis=0) * 100
ct_per_100c["n_cells"] = n_cells
print("\n=== Reads per 100 cells per cluster (cell-density-normalized) ===")
print(ct_per_100c.round(2).to_string())
ct_per_100c.to_csv(BASE / "novel_isoforms_reads_per_100_cells.tsv", sep="\t")

print("\nWrote:")
for f in ["novel_isoforms_per_cluster_pct.tsv",
          "novel_isoforms_per_cluster_counts.tsv",
          "novel_isoforms_reads_per_100_cells.tsv"]:
    print(f"  {BASE / f}")
