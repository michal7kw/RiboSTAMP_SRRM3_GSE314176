#!/usr/bin/env python3
"""Count cells per cluster in the long-read dataset.

Reports two distinct numbers:
  1) cells in the published GSE314176 long-read whitelist (per mouse + pooled),
     for each Cell Assignments Grouped (high-level: Astro, Oligo, DG, CA1, CA3,
     GABA, Endo, ...) and each Cell Assignment (subtype, e.g. Neuron - DG_5).
  2) cells we actually had ≥1 read at the Srrm3 locus AND a matched cell
     barcode in step 03 (these are the cells contributing to per-cluster PSI).
"""
from collections import Counter
from pathlib import Path

import pandas as pd

META_DIR = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata")
PER_READ = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/per_read/all_samples.per_read.tsv")

SRR_TO_GSM = {
    "SRR36480452": "GSM9380801",   # mouse3
    "SRR36480453": "GSM9380800",   # mouse2
    "SRR36480454": "GSM9380799",   # mouse1
}

OBS_FILES = {
    "GSM9380799 (mouse1)": META_DIR / "GSM9380799_longread_normed_counts_transcript_adata_mouse1__obs.tsv",
    "GSM9380800 (mouse2)": META_DIR / "GSM9380800_longread_normed_counts_transcript_adata_mouse2__obs.tsv",
    "GSM9380801 (mouse3)": META_DIR / "GSM9380801_longread_normed_counts_transcript_adata_mouse3__obs.tsv",
}

# ---------- 1. published whitelist counts ----------
print("=" * 70)
print("1. CELLS IN THE PUBLISHED LONG-READ WHITELIST")
print("=" * 70)

pooled_grouped = Counter()
pooled_assignment = Counter()

for label, path in OBS_FILES.items():
    df = pd.read_csv(path, sep="\t")
    grouped = df["Cell Assignments Grouped"].value_counts()
    print(f"\n>>> {label}  (n={len(df)} cells)")
    for ct, n in grouped.items():
        print(f"   {ct:<30s}{n:>6d}")
    pooled_grouped.update(df["Cell Assignments Grouped"].dropna())
    pooled_assignment.update(df["Cell Assignment"].dropna())

print("\n--- POOLED (3 mice combined) — Cell Assignments Grouped ---")
total = sum(pooled_grouped.values())
for ct, n in sorted(pooled_grouped.items(), key=lambda x: -x[1]):
    print(f"   {ct:<30s}{n:>6d}   ({100*n/total:5.1f}%)")
print(f"   {'TOTAL':<30s}{total:>6d}")

print("\n--- POOLED — Cell Assignment (subtypes) ---")
for ct, n in sorted(pooled_assignment.items(), key=lambda x: -x[1]):
    print(f"   {ct:<30s}{n:>6d}")

# ---------- 2. cells detected at Srrm3 ----------
print("\n" + "=" * 70)
print("2. CELLS WITH ≥1 READ AT Srrm3 + MATCHED CELL BARCODE")
print("   (i.e. cells contributing to per-cluster PSI in step 04)")
print("=" * 70)

per_read = pd.read_csv(PER_READ, sep="\t")
print(f"\nTotal per-read rows: {len(per_read):,}")
print(f"Rows with matched CB:  {per_read['cb_matched'].notna().sum():,}")

# bare barcode = strip _sample suffix
def bare(bc):
    if pd.isna(bc):
        return None
    return str(bc).split("_", 1)[0]

# Build (mouse, bare_barcode) -> cluster map from obs files
cluster_map = {}
for label, path in OBS_FILES.items():
    gsm = label.split()[0]
    df = pd.read_csv(path, sep="\t")
    for _, row in df.iterrows():
        bc = bare(row["barcode"])
        cluster_map[(gsm, bc)] = (row["Cell Assignment"], row["Cell Assignments Grouped"])

# Annotate per_read with cluster
per_read["mouse_gsm"] = per_read["sample"].map(SRR_TO_GSM)
per_read["bare_cb"] = per_read["cb_matched"].apply(bare)

detected_cells = (
    per_read.dropna(subset=["bare_cb"])
            .drop_duplicates(["mouse_gsm", "bare_cb"])
)
print(f"\nUnique cells detected at Srrm3 (matched CB): {len(detected_cells):,}")

# Map to clusters
clusters_grouped = []
clusters_subtype = []
unmapped = 0
for _, row in detected_cells.iterrows():
    key = (row["mouse_gsm"], row["bare_cb"])
    if key in cluster_map:
        sub, grp = cluster_map[key]
        clusters_subtype.append(sub)
        clusters_grouped.append(grp)
    else:
        unmapped += 1

print(f"Cells with cluster label:    {len(clusters_grouped):,}")
print(f"Cells without cluster label: {unmapped:,}  (matched CB but not in published whitelist)")

print("\n--- Detected cells, by Cell Assignments Grouped (POOLED) ---")
g = Counter(clusters_grouped)
total_g = sum(g.values())
for ct, n in sorted(g.items(), key=lambda x: -x[1]):
    print(f"   {ct:<30s}{n:>6d}   ({100*n/total_g:5.1f}%)")
print(f"   {'TOTAL':<30s}{total_g:>6d}")

print("\n--- Detected cells, by Cell Assignments Grouped — per mouse ---")
for label, path in OBS_FILES.items():
    gsm = label.split()[0]
    sub = detected_cells[detected_cells["mouse_gsm"] == gsm]
    counts = Counter()
    for _, row in sub.iterrows():
        key = (gsm, row["bare_cb"])
        if key in cluster_map:
            counts[cluster_map[key][1]] += 1
    print(f"\n>>> {label}")
    tot = sum(counts.values())
    for ct, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"   {ct:<30s}{n:>6d}")
    print(f"   {'subtotal':<30s}{tot:>6d}")

print("\n--- Detected cells, by Cell Assignment subtype (POOLED, top 25) ---")
s = Counter(clusters_subtype)
for ct, n in sorted(s.items(), key=lambda x: -x[1])[:25]:
    print(f"   {ct:<30s}{n:>6d}")

# ---------- 3. cluster capture rate (detected / whitelist) ----------
print("\n" + "=" * 70)
print("3. PER-CLUSTER CAPTURE RATE — (detected at Srrm3) / (whitelist)")
print("=" * 70)
print(f"   {'Cluster':<28s}{'whitelist':>12s}{'detected':>12s}{'rate':>8s}")
for ct, n_wl in sorted(pooled_grouped.items(), key=lambda x: -x[1]):
    n_det = g.get(ct, 0)
    rate = 100 * n_det / n_wl if n_wl else 0
    print(f"   {ct:<28s}{n_wl:>12d}{n_det:>12d}{rate:>7.1f}%")
