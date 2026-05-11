#!/usr/bin/env python3
"""Surface B (cheap variant) — re-classify long-read reads against the
short-read whitelist so we can recover OPC + DG-Immature progenitors.

Why this works
--------------
10x Chromium cell barcodes are universal: the same physical droplet has
the same CB sequence in both long-read and short-read libraries.

The published GSE314176 LONG-READ whitelist filtered out OPC + DG-Immature
during their QC (those cells had too few transcripts for confident isoform
assignment). The published SHORT-READ whitelist DOES contain those cells.

Plan:
  1. Load short-read obs.tsv per mouse → barcode → cell_type lookup
     (filtered to OPC + DG-Immature)
  2. Load long-read per_read.tsv (already has cb_extracted column from
     polyA-anchored extraction)
  3. For each per_read row, look up cb_extracted against the short-read
     OPC/DG-Immature barcode set (allowing Hamming-1 fallback to absorb
     PacBio sequencing errors in the CB region)
  4. If matched, count INCLUSION/SKIP/UNINFORMATIVE per progenitor cell type
  5. Compute Wilson 95% CI per group and emit verdict markdown

What this DOES vs what Surface B (HPC) would do
-----------------------------------------------
This works on the LONG-READ BAMs we already have aligned. If the long-
read pipeline filtered OPC/DG cells out at the read level (not just at
the cell-whitelist level), this analysis will find zero matches and we
need the full Surface B (cellranger on raw 10x FASTQ).

Most likely: OPC/DG cells were filtered at the cell-whitelist level
(not at read level). So their reads ARE in the long-read BAMs but had
no cluster label assigned. This script recovers them.

Outputs:
  data/progenitor_psi/per_progenitor_cell_psi.tsv  — per cell × class
  data/progenitor_psi/per_cluster_psi.tsv          — per progenitor cluster
  data/progenitor_psi/progenitor_verdict.md        — written verdict
"""
from __future__ import annotations
import os
import sys
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd
from scipy.stats import binomtest

BASE_DIR = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data")
META_DIR = BASE_DIR / "metadata"
PER_READ = BASE_DIR / "targeted_psi" / "per_read" / "all_samples.per_read.tsv"
OUT_DIR = BASE_DIR / "progenitor_psi"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# SRR (long-read) → mouse number, since short-read obs is keyed by mouse
SRR_TO_MOUSE = {
    "SRR36480452": "mouse3",   # long-read = GSM9380801
    "SRR36480453": "mouse2",   # long-read = GSM9380800
    "SRR36480454": "mouse1",   # long-read = GSM9380799
}
SHORTREAD_GSM_PER_MOUSE = {
    "mouse1": "GSM9380796",
    "mouse2": "GSM9380797",
    "mouse3": "GSM9380798",
}

# Cell types of interest — the progenitor populations
PROGENITOR_TYPES = {"Glia - OPC", "Neuron - DG - Immature"}

# Hamming-1 tolerance lookup (sequencing errors in CB)
def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def build_hamming1_index(barcodes):
    """For Hamming-1 lookup: build a dict of all 1-edit neighbors → original."""
    bc_set = set(barcodes)
    h1 = dict()
    for b in barcodes:
        for i in range(len(b)):
            for nt in "ACGT":
                if nt == b[i]:
                    continue
                neighbor = b[:i] + nt + b[i+1:]
                # If neighbor is itself a real barcode, exact match wins;
                # otherwise this is a 1-edit lookup target
                if neighbor not in bc_set:
                    if neighbor in h1 and h1[neighbor] != b:
                        # Ambiguous — drop it (don't risk wrong assignment)
                        h1[neighbor] = None
                    else:
                        h1[neighbor] = b
    return h1, bc_set

# -----------------------------------------------------------------------------
# Step 1: Load progenitor barcodes per mouse from short-read obs.tsv
# -----------------------------------------------------------------------------
print("=" * 70)
print("Step 1: Load progenitor cells from short-read whitelists")
print("=" * 70)

progenitor_bc_per_mouse = {}  # mouse → set of barcodes
progenitor_bc_to_celltype = {}  # (mouse, barcode) → cell_type

for mouse, gsm in SHORTREAD_GSM_PER_MOUSE.items():
    obs_path = list(META_DIR.glob(f"{gsm}_shortread_*_obs.tsv"))[0]
    df = pd.read_csv(obs_path, sep="\t")

    # Pull the cell type column
    if "Cell Assignments Grouped" not in df.columns:
        sys.exit(f"ERROR: no 'Cell Assignments Grouped' column in {obs_path}")

    progenitors = df[df["Cell Assignments Grouped"].isin(PROGENITOR_TYPES)].copy()
    # Strip ANY suffix: long-read uses "_sampleN", short-read uses "-1".
    # Keep only leading ACGT characters — barcodes are pure ACGT.
    _BC_RE = re.compile(r"^([ACGT]+)")
    def _strip(s):
        m = _BC_RE.match(str(s))
        return m.group(1) if m else str(s)
    progenitors["bare_barcode"] = progenitors["barcode"].astype(str).map(_strip)

    print(f"\n  {mouse} ({gsm}):")
    n_total = len(df)
    n_progenitors = len(progenitors)
    by_type = progenitors["Cell Assignments Grouped"].value_counts()
    print(f"    Total cells:               {n_total:,}")
    print(f"    Progenitor cells found:    {n_progenitors:,}")
    for ct, n in by_type.items():
        print(f"      {ct}: {n}")

    bcs = set(progenitors["bare_barcode"])
    progenitor_bc_per_mouse[mouse] = bcs
    for _, row in progenitors.iterrows():
        progenitor_bc_to_celltype[(mouse, row["bare_barcode"])] = row["Cell Assignments Grouped"]

# Pool across mice
all_progenitor_bcs = set()
for s in progenitor_bc_per_mouse.values():
    all_progenitor_bcs |= s
print(f"\nTotal unique progenitor barcodes across 3 mice: {len(all_progenitor_bcs)}")

# -----------------------------------------------------------------------------
# Step 2: Build Hamming-1 lookup index per mouse
# -----------------------------------------------------------------------------
print("\nStep 2: Build Hamming-1 lookup indices")
hamming1_index_per_mouse = {}
for mouse, bcs in progenitor_bc_per_mouse.items():
    h1, exact = build_hamming1_index(bcs)
    hamming1_index_per_mouse[mouse] = (exact, h1)
    print(f"  {mouse}: {len(exact)} exact + {len(h1)} hamming-1 candidates")

# -----------------------------------------------------------------------------
# Step 3: Stream long-read per_read.tsv and lookup cb_extracted
# -----------------------------------------------------------------------------
print("\nStep 3: Stream long-read per_read.tsv and lookup matches")

if not PER_READ.exists():
    sys.exit(f"ERROR: {PER_READ} not found")

# Per-cell counts (mouse, bare_barcode) → {classification → count}
per_cell_counts = defaultdict(lambda: defaultdict(int))

# Aggregate counts (mouse, cell_type) → {classification → count}
per_cluster_counts = defaultdict(lambda: defaultdict(int))

# Track everything for diagnostics
n_total_rows = 0
n_with_cb_extracted = 0
n_exact_match = 0
n_hamming1_match = 0
n_unmatched = 0
match_orientation = defaultdict(int)

for chunk in pd.read_csv(PER_READ, sep="\t", chunksize=200_000):
    for _, row in chunk.iterrows():
        n_total_rows += 1
        srr = row["sample"]
        mouse = SRR_TO_MOUSE.get(srr)
        if mouse is None:
            continue
        cb = row["cb_extracted"]
        if pd.isna(cb) or cb == "" or cb == "NA":
            continue
        n_with_cb_extracted += 1

        exact, h1 = hamming1_index_per_mouse[mouse]
        matched = None
        if cb in exact:
            matched = cb
            n_exact_match += 1
        elif cb in h1:
            corrected = h1[cb]
            if corrected is not None:
                matched = corrected
                n_hamming1_match += 1
        if matched is None:
            n_unmatched += 1
            continue

        cell_type = progenitor_bc_to_celltype.get((mouse, matched))
        if cell_type is None:
            continue  # shouldn't happen

        cls = row["classification"]
        per_cell_counts[(mouse, matched, cell_type)][cls] += 1
        per_cluster_counts[(mouse, cell_type)][cls] += 1
        match_orientation[row.get("cb_orientation", "?")] += 1

print(f"\n  per_read rows scanned:                    {n_total_rows:>10,}")
print(f"  rows with cb_extracted (non-empty):       {n_with_cb_extracted:>10,}")
print(f"  exact match to progenitor barcode:        {n_exact_match:>10,}")
print(f"  hamming-1 match:                          {n_hamming1_match:>10,}")
print(f"  unmatched:                                {n_unmatched:>10,}")
print(f"  total matched:                            {n_exact_match + n_hamming1_match:>10,}")

if n_exact_match + n_hamming1_match == 0:
    print("\nNO MATCHES FOUND — likely the OPC/DG-Immature cells' reads were filtered")
    print("at the read level by the long-read pipeline, not just at the cell-whitelist")
    print("level. To proceed: run the full Surface B pipeline (cellranger on raw 10x FASTQ).")
    sys.exit(0)

# -----------------------------------------------------------------------------
# Step 4: Aggregate per cluster (across mice) + compute PSI + Wilson CIs
# -----------------------------------------------------------------------------
print("\nStep 4: Aggregate per cluster + compute PSI")

# First: pool per-mouse counts to per-cluster across mice
pooled_per_cluster = defaultdict(lambda: defaultdict(int))
n_cells_per_cluster = defaultdict(int)
for (mouse, ct), counts in per_cluster_counts.items():
    for cls, n in counts.items():
        pooled_per_cluster[ct][cls] += n
# Count unique cells (matched barcodes) per cluster
unique_cells_per_cluster = defaultdict(set)
for (mouse, bc, ct), _ in per_cell_counts.items():
    unique_cells_per_cluster[ct].add((mouse, bc))
for ct, cells in unique_cells_per_cluster.items():
    n_cells_per_cluster[ct] = len(cells)

def wilson(num, den):
    if den == 0:
        return None, None, None
    res = binomtest(int(num), int(den))
    lo, hi = res.proportion_ci(method="wilson")
    return num/den, lo, hi

cluster_rows = []
for ct in sorted(pooled_per_cluster.keys()):
    counts = pooled_per_cluster[ct]
    inc = counts.get("INCLUSION", 0)
    skip = counts.get("SKIP", 0)
    uninform = counts.get("UNINFORMATIVE", 0)
    total_inform = inc + skip
    psi, lo, hi = wilson(inc, total_inform)
    cluster_rows.append({
        "cell_type": ct,
        "n_cells": n_cells_per_cluster[ct],
        "n_inclusion": inc,
        "n_skip": skip,
        "n_uninformative": uninform,
        "n_informative": total_inform,
        "PSI": psi,
        "PSI_95_CI_low": lo,
        "PSI_95_CI_high": hi,
    })

cluster_df = pd.DataFrame(cluster_rows)

# Add a "all progenitors pooled" row
all_inc = sum(r["n_inclusion"] for r in cluster_rows)
all_skip = sum(r["n_skip"] for r in cluster_rows)
all_uninform = sum(r["n_uninformative"] for r in cluster_rows)
all_inform = all_inc + all_skip
psi_all, lo_all, hi_all = wilson(all_inc, all_inform)
all_n_cells = sum(r["n_cells"] for r in cluster_rows)
cluster_df = pd.concat([cluster_df, pd.DataFrame([{
    "cell_type": "ALL_PROGENITORS_POOLED",
    "n_cells": all_n_cells,
    "n_inclusion": all_inc,
    "n_skip": all_skip,
    "n_uninformative": all_uninform,
    "n_informative": all_inform,
    "PSI": psi_all,
    "PSI_95_CI_low": lo_all,
    "PSI_95_CI_high": hi_all,
}])], ignore_index=True)

print()
print(cluster_df.to_string(index=False))

# Save outputs
PER_CLUSTER_TSV = OUT_DIR / "per_cluster_psi.tsv"
cluster_df.to_csv(PER_CLUSTER_TSV, sep="\t", index=False)
print(f"\nWrote: {PER_CLUSTER_TSV}")

# Per-cell table
per_cell_rows = []
for (mouse, bc, ct), counts in per_cell_counts.items():
    per_cell_rows.append({
        "mouse": mouse,
        "barcode": bc,
        "cell_type": ct,
        "n_inclusion": counts.get("INCLUSION", 0),
        "n_skip": counts.get("SKIP", 0),
        "n_uninformative": counts.get("UNINFORMATIVE", 0),
    })
per_cell_df = pd.DataFrame(per_cell_rows)
PER_CELL_TSV = OUT_DIR / "per_progenitor_cell_psi.tsv"
per_cell_df.to_csv(PER_CELL_TSV, sep="\t", index=False)
print(f"Wrote: {PER_CELL_TSV}")

# -----------------------------------------------------------------------------
# Step 5: Verdict markdown
# -----------------------------------------------------------------------------
print("\nStep 5: Verdict")

verdict_lines = []
verdict_lines.append("# Surface B (cheap variant) — progenitor PSI verdict\n")
verdict_lines.append(f"**Date:** 2026-04-29\n")
verdict_lines.append(f"**Method:** cross-library cell-barcode matching (long-read reads → short-read whitelist for OPC + DG-Immature)\n")
verdict_lines.append(f"**Source:** [`pc_pipeline/targeted_psi/07_progenitor_psi_via_shortread_cb.py`](../../pc_pipeline/targeted_psi/07_progenitor_psi_via_shortread_cb.py)\n\n")

verdict_lines.append("## Match statistics\n\n")
verdict_lines.append(f"| Metric | Count |\n|---|---:|\n")
verdict_lines.append(f"| per_read rows scanned | {n_total_rows:,} |\n")
verdict_lines.append(f"| With non-empty `cb_extracted` | {n_with_cb_extracted:,} |\n")
verdict_lines.append(f"| Exact match to progenitor barcode | {n_exact_match:,} |\n")
verdict_lines.append(f"| Hamming-1 match | {n_hamming1_match:,} |\n")
verdict_lines.append(f"| Total matched to progenitor | {n_exact_match + n_hamming1_match:,} |\n\n")

verdict_lines.append("## Per-cluster PSI\n\n")
verdict_lines.append(cluster_df.to_markdown(index=False) + "\n\n")

# Compare to mature types
mature_max_psi_pct = 0.05  # from per_cluster_psi.tsv Wilson upper CI
verdict_lines.append("## Verdict\n\n")
if all_inform == 0:
    verdict = "INCONCLUSIVE"
    interp = (
        "ZERO informative reads from progenitor barcodes were found in the long-read library. "
        "This means OPC + DG-Immature cells' reads were filtered at the read level by the "
        "long-read pipeline (not just at the cell-whitelist level). To test the progenitor "
        "hypothesis, the full Surface B pipeline (cellranger on raw 10x FASTQ) is required."
    )
elif psi_all is not None and psi_all > 0.05:
    verdict = "PROGENITOR HYPOTHESIS CONFIRMED"
    interp = (
        f"Pooled progenitor PSI = {psi_all:.4f} ({psi_all*100:.2f}%) — substantially "
        f"higher than mature hippocampal types (all <0.05%). The cassette IS enriched "
        f"in progenitor populations of mature hippocampus."
    )
elif psi_all is not None and psi_all > 0:
    verdict = "INTERMEDIATE"
    interp = (
        f"Pooled progenitor PSI = {psi_all:.4f} ({psi_all*100:.4f}%) — non-zero but "
        f"not dramatically enriched relative to mature cell types. Underpowered or "
        f"genuinely modest progenitor enrichment."
    )
else:
    verdict = "PROGENITOR HYPOTHESIS REJECTED"
    interp = (
        f"Pooled progenitor PSI = 0% (Wilson upper CI {hi_all*100:.4f}%) — the cassette "
        f"is NOT enriched in OPC or DG-Immature cells either. The cassette is a cell-line "
        f"/ proliferation-restricted phenomenon with no detectable in vivo enrichment in "
        f"the captured progenitor populations."
    )
verdict_lines.append(f"**{verdict}**\n\n{interp}\n\n")

verdict_lines.append(f"## Comparison to other contexts\n\n")
verdict_lines.append("| Context | PSI | Source |\n|---|---:|---|\n")
verdict_lines.append("| Proliferating N2A — Parental (lab anchor) | 57% | `90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md` |\n")
verdict_lines.append("| Proliferating N2A — Parental (our calc) | 50.6% | `data/lab_bam_validation/lab_n2a_psi_per_group.tsv` |\n")
verdict_lines.append("| Proliferating N2A — KO | 5% / 2.7% | same |\n")
verdict_lines.append("| Differentiated cell-line neurons (NT/B15/B60) | 0% | `data/sr_junction_psi/results/cell_line_validation_verdict_final.md` |\n")
verdict_lines.append("| Adult hippocampus mature types (7) | <0.05% | `data/targeted_psi/results/per_cluster_psi.tsv` |\n")
psi_pct_str = "—" if psi_all is None else f"{psi_all*100:.4f}%"
verdict_lines.append(f"| **OPC + DG-Immature progenitors (this analysis)** | **{psi_pct_str}** | `data/progenitor_psi/per_cluster_psi.tsv` |\n\n")

VERDICT_MD = OUT_DIR / "progenitor_verdict.md"
with open(VERDICT_MD, "w") as f:
    f.writelines(verdict_lines)
print(f"\nWrote: {VERDICT_MD}")
print(f"\nFinal verdict: {verdict}")
print(interp)
