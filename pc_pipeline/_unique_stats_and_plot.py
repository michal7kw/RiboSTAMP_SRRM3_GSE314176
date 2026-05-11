#!/usr/bin/env python3
"""Stats + plot on unique-only counts (rigorous baseline)."""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

BASE = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/results")
counts = pd.read_csv(BASE / "novel_isoforms_per_cluster_UNIQUE_counts.tsv", sep="\t", index_col=0)
pct = pd.read_csv(BASE / "novel_isoforms_per_cluster_UNIQUE_pct.tsv", sep="\t", index_col=0)

order = ["BBB - Endo", "Glia - Astro", "Glia - Oligo",
         "Neuron - DG", "Neuron - CA1", "Neuron - CA3", "Neuron - GABA"]
order = [o for o in order if o in counts.index]
counts = counts.reindex(order)
pct = pct.reindex(order)

# Stats
non_neuron = ["BBB - Endo", "Glia - Astro", "Glia - Oligo"]
neuron = ["Neuron - DG", "Neuron - CA1", "Neuron - CA3", "Neuron - GABA"]
novel_cols = [c for c in counts.columns if c.startswith("transcript")]
canon = "ENSMUST00000005077.6"

n_novel_neu = counts.loc[neuron, novel_cols].values.sum()
n_canon_neu = counts.loc[neuron, canon].sum()
n_novel_glia = counts.loc[non_neuron, novel_cols].values.sum()
n_canon_glia = counts.loc[non_neuron, canon].sum()

table = pd.DataFrame([[n_novel_neu, n_canon_neu], [n_novel_glia, n_canon_glia]],
                     index=["neuron", "non-neuron"],
                     columns=["novel (unique)", "canonical Srrm3-201 (unique)"])
print("=== Chi-square (unique-only counts) ===")
print(table)
chi2, p, dof, exp = chi2_contingency(table.values)
print(f"chi2 = {chi2:.2f}, dof = {dof}, p = {p:.3e}")
or_ = (n_novel_neu * n_canon_glia) / (n_novel_glia * n_canon_neu)
print(f"odds ratio = {or_:.2f}")

print("\n=== Per-novel-transcript Fisher (unique-only) ===")
for nv in novel_cols:
    a = counts.loc[neuron, nv].sum()
    b = counts.loc[neuron, canon].sum()
    c = counts.loc[non_neuron, nv].sum()
    d = counts.loc[non_neuron, canon].sum()
    or2, p2 = fisher_exact([[a, b], [c, d]], alternative="greater")
    print(f"  {nv}: neuron={a}/{b}, non-neuron={c}/{d}, OR={or2:.2f}, p={p2:.3e}")

# Stacked bar (zoomed in on the novel fraction; canonical bar always >90%)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5),
                                gridspec_kw={"width_ratios": [1, 1.5]})

iso_order = ["ENSMUST00000005077.6", "ENSMUST00000144211.1",
             "transcript185.chr5.nnic", "transcript11.chr5.nnic",
             "transcript157.chr5.nnic"]
iso_order = [c for c in iso_order if c in pct.columns]
labels = {
    "ENSMUST00000005077.6":    "Srrm3-201 (canonical)",
    "ENSMUST00000144211.1":    "Srrm3-204 (ref)",
    "transcript185.chr5.nnic": "novel t185 (skips Ex15)",
    "transcript11.chr5.nnic":  "novel t11 (NMD)",
    "transcript157.chr5.nnic": "novel t157 (extra intron)",
}
colors = {
    "ENSMUST00000005077.6":    "#5A5A5A",
    "ENSMUST00000144211.1":    "#1D7874",
    "transcript185.chr5.nnic": "#E63946",
    "transcript11.chr5.nnic":  "#457B9D",
    "transcript157.chr5.nnic": "#D08C34",
}

# Left: full stacked bar (canonical dominates)
bottom = np.zeros(len(pct))
for c in iso_order:
    ax1.bar(pct.index, pct[c], bottom=bottom, color=colors[c],
            label=labels[c], edgecolor="white", linewidth=0.5)
    bottom += pct[c].values
ax1.set_ylabel("% of cluster Srrm3 reads (unique only)")
ax1.set_title("Full view (0–100%)")
ax1.set_ylim(0, 100)
plt.setp(ax1.get_xticklabels(), rotation=20, ha="right")

# Right: zoomed in on the novel-fraction (0-8%)
novel_only = [c for c in iso_order if c.startswith("transcript")] + ["ENSMUST00000144211.1"]
bottom = np.zeros(len(pct))
for c in novel_only:
    ax2.bar(pct.index, pct[c], bottom=bottom, color=colors[c],
            label=labels[c], edgecolor="white", linewidth=0.5)
    bottom += pct[c].values
ax2.set_ylabel("% of cluster Srrm3 reads (novel + Srrm3-204 only)")
ax2.set_title("Zoomed: novel + ref-Srrm3-204 fraction\n"
              "(canonical Srrm3-201 omitted)")
ax2.set_ylim(0, max(8, bottom.max() * 1.1))
plt.setp(ax2.get_xticklabels(), rotation=20, ha="right")
ax2.legend(loc="upper left", fontsize=9)

plt.suptitle("Srrm3 isoform composition per hippocampal cell type — UNIQUE assignments only\n"
             f"(P25 mouse, GSE314176, n={int(counts['TOTAL_unique'].sum()):,} unique-assigned reads)",
             y=1.02)
plt.tight_layout()
out = BASE / "novel_isoforms_per_cluster_UNIQUE.png"
plt.savefig(out, dpi=180, bbox_inches="tight")
print(f"\nWrote {out}")
