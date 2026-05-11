#!/usr/bin/env python3
"""
Per-cluster Srrm3 isoform composition: stacked barplot + chi-square test of
neuronal-vs-non-neuronal enrichment for the novel isoforms.
"""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

BASE = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/results")

counts = pd.read_csv(BASE / "novel_isoforms_per_cluster_counts.tsv", sep="\t", index_col=0)
pct = pd.read_csv(BASE / "novel_isoforms_per_cluster_pct.tsv", sep="\t", index_col=0)

# Reorder cell types: glia first, neurons grouped, endo last
order = [
    "BBB - Endo", "Glia - Astro", "Glia - Oligo",
    "Neuron - DG", "Neuron - CA1", "Neuron - CA3", "Neuron - GABA",
]
order = [o for o in order if o in counts.index]
counts = counts.reindex(order)
pct = pct.reindex(order)

# --- Stacked bar of % composition ---
iso_cols = [c for c in [
    "ENSMUST00000005077.6",
    "ENSMUST00000144211.1",
    "transcript185.chr5.nnic",
    "transcript11.chr5.nnic",
    "transcript157.chr5.nnic",
] if c in pct.columns]
labels = {
    "ENSMUST00000005077.6":    "Srrm3-201 (canonical)",
    "ENSMUST00000144211.1":    "Srrm3-204 (ref)",
    "transcript185.chr5.nnic": "novel t185 (skips Ex15)",
    "transcript11.chr5.nnic":  "novel t11 (11 exons)",
    "transcript157.chr5.nnic": "novel t157 (extra intron)",
}
colors = {
    "ENSMUST00000005077.6":    "#5A5A5A",   # gray
    "ENSMUST00000144211.1":    "#1D7874",   # teal
    "transcript185.chr5.nnic": "#E63946",   # warm red
    "transcript11.chr5.nnic":  "#457B9D",   # blue
    "transcript157.chr5.nnic": "#D08C34",   # ochre
}

fig, ax = plt.subplots(figsize=(8, 5.5))
bottom = np.zeros(len(pct))
for c in iso_cols:
    ax.bar(pct.index, pct[c], bottom=bottom,
           color=colors[c], label=labels[c], edgecolor="white", linewidth=0.5)
    bottom += pct[c].values
ax.set_ylabel("% of cluster's Srrm3 reads")
ax.set_title("Srrm3 isoform composition per hippocampal cell type\n"
             "(P25 mouse, GSE314176 long-read, n=24,067 reads)")
ax.set_ylim(0, 100)
ax.legend(loc="upper right", bbox_to_anchor=(1.42, 1.0), fontsize=9)
plt.xticks(rotation=20, ha="right")
plt.tight_layout()
out_png = BASE / "novel_isoforms_per_cluster.png"
plt.savefig(out_png, dpi=180, bbox_inches="tight")
plt.close()
print(f"Wrote {out_png}")

# --- Statistical test: novel transcripts enriched in neurons? ---
non_neuron = ["BBB - Endo", "Glia - Astro", "Glia - Oligo"]
neuron     = ["Neuron - DG", "Neuron - CA1", "Neuron - CA3", "Neuron - GABA"]

print("\n=== Chi-square test: novel-vs-canonical × neuron-vs-non-neuron ===")
print("(novel = transcript11 + transcript157 + transcript185)")

novel_cols = [c for c in iso_cols if c.startswith("transcript")]
canonical_cols = [c for c in iso_cols if not c.startswith("transcript")]

n_novel_neu      = counts.loc[neuron, novel_cols].values.sum()
n_canonical_neu  = counts.loc[neuron, canonical_cols].values.sum()
n_novel_glia     = counts.loc[non_neuron, novel_cols].values.sum()
n_canonical_glia = counts.loc[non_neuron, canonical_cols].values.sum()

table = pd.DataFrame(
    [[n_novel_neu, n_canonical_neu], [n_novel_glia, n_canonical_glia]],
    index=["neuron", "non-neuron"], columns=["novel", "canonical/known-ref"]
)
print(table)
chi2, p, dof, exp = chi2_contingency(table.values)
print(f"\nchi2 = {chi2:.2f},  dof = {dof},  p = {p:.3e}")
odds_ratio = (n_novel_neu * n_canonical_glia) / (n_novel_glia * n_canonical_neu)
print(f"odds ratio (novel|neuron vs novel|non-neuron) = {odds_ratio:.2f}")

# Per-novel-transcript Fisher
print("\n=== Per-novel-transcript Fisher exact (each novel vs canonical Srrm3-201) ===")
canon = "ENSMUST00000005077.6"
for nv in novel_cols:
    a = counts.loc[neuron, nv].sum()
    b = counts.loc[neuron, canon].sum()
    c = counts.loc[non_neuron, nv].sum()
    d = counts.loc[non_neuron, canon].sum()
    or_, p2 = fisher_exact([[a, b], [c, d]], alternative="greater")
    print(f"  {labels[nv]}: neuron-novel/canon={a}/{b}, non-neuron-novel/canon={c}/{d},  OR={or_:.2f},  p={p2:.3e}")

# --- Also save markdown summary ---
md = BASE / "ISOFORM_COMPOSITION.md"
with open(md, "w") as fh:
    fh.write("# Srrm3 isoform composition per cell type\n\n")
    fh.write("## Cluster-normalized % composition\n\n")
    fh.write("```\n")
    cols_pretty = pct.copy()
    cols_pretty.columns = [labels.get(c, c) if c in labels else c for c in cols_pretty.columns]
    fh.write(cols_pretty.round(2).to_string() + "\n```\n\n")
    fh.write(f"## Test: novel isoforms enriched in neurons?\n\n")
    fh.write("```\n")
    fh.write(table.to_string() + "\n")
    fh.write(f"\nchi2 = {chi2:.2f},  dof = {dof},  p = {p:.3e}\n")
    fh.write(f"odds ratio (novel|neuron) = {odds_ratio:.2f}\n```\n\n")
    fh.write("Interpretation: novel Srrm3 isoforms collectively account for ")
    fh.write(f"~{table['novel']['neuron']/table.loc['neuron'].sum()*100:.0f}% of Srrm3 reads in hippocampal neurons ")
    fh.write(f"vs ~{table['novel']['non-neuron']/table.loc['non-neuron'].sum()*100:.0f}% in glia/endothelial.\n")
print(f"\nWrote {md}")
