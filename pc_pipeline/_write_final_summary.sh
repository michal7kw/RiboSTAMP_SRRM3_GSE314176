#!/bin/bash
set -euo pipefail
RES=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/results
SUM="$RES/SUMMARY.md"

cat > "$SUM" <<'EOF'
# Path A — De-novo Srrm3 isoform discovery (GSE314176 hippocampus)

**Date:** 2026-05-08
**Pipeline:** `pc_pipeline` Path A (local-PC variant) — IsoQuant 3.12.2 (`sensitive_pacbio`) + SQANTI3 v5.4
**Compute:** ~25 min total local (i5-14600K, 16 threads, WSL2)
**Full PI-facing writeup:** `../../../docs/PI_RESPONSE_04_novel_srrm3_isoforms.md`
**Methodology rationale:** `../../../docs/LEARNING_11_isoform_discovery_fastpath.md`
**Reproducible commands:** `../../../docs/NOVEL_ISOFORMS_RUNBOOK.md`

## TL;DR

> **Three high-confidence novel Srrm3 isoforms discovered in P25 mouse hippocampus.** All three carry at least one novel splice junction (NNC), pass intra-priming filters (perc_A_downstream < 30%), and replicate across all 3 mice. **All three are 4-5× enriched in neurons relative to glia/endothelial cells (chi² = 192, p ≈ 10⁻⁴³, OR = 4.62, on rigorous unique-only assignments).** Two of them (`t11`, `t157`) encode the **same 95 aa truncated protein** through different splice paths (one predicted NMD); the third (`t185`) is **non-coding** and spans the lab's Ex15 cassette region but skips it via a 3.8-kb intron.

## The 3 high-confidence novel isoforms

| Transcript | Cat. | Exons | Length | Novel feature | Coding | ORF | NMD |
|---|---|---:|---:|---|---|---:|---|
| **transcript185.chr5.nnic** | NNC | 4 | 2802 | spans Ex15 region but skips it; 4-exon 3'-end form | **non-coding** | — | — |
| **transcript11.chr5.nnic** | NNC | 11 | 1974 | terminal exon shift; 11-exon long form | coding | 95 aa | **TRUE** |
| **transcript157.chr5.nnic** | NNC | 4 | 742 | extra novel intron (AT-AC, possible U12) | coding | 95 aa (same as t11) | FALSE |

## Cluster-resolved expression — UNIQUE-only counts (rigorous baseline)

```
                Srrm3-201   Srrm3-204   t185 novel   t11 novel   t157 novel
BBB - Endo       99.6%       0.0%        0.1%         0.2%        0.1%
Glia - Astro     99.1%       0.1%        0.2%         0.5%        0.1%
Glia - Oligo     96.1%       0.5%        1.1%         1.9%        0.5%
Neuron - DG      94.7%       0.4%        0.9%         3.1%        0.9%
Neuron - CA1     93.2%       0.5%        1.9%         3.5%        0.9%
Neuron - CA3     95.7%       0.4%        0.2%         2.6%        1.0%
Neuron - GABA    94.1%       0.2%        1.7%         3.6%        0.5%
```

Each novel isoform individually accounts for ~1-4% of cluster-level Srrm3 reads in neurons, dominated by canonical Srrm3-201 (>93%). The neuron enrichment is qualitative (frequency in neurons / frequency in non-neurons), not absolute dominance.

## Statistical test — neuron enrichment (UNIQUE-only)

| Group | Novel reads | Canonical Srrm3-201 reads |
|---|---:|---:|
| **Neurons** (CA1, CA3, DG, GABA) | 152 | 2,650 |
| **Non-neurons** (Endo, Astro, Oligo) | 136 | 10,962 |

**chi² = 192.36 (dof=1), p = 9.7 × 10⁻⁴⁴**, **odds ratio = 4.62**

Per novel transcript (Fisher exact, neuron vs non-neuron, vs canonical Srrm3-201):

| Transcript | OR | p (one-sided) |
|---|---:|---:|
| transcript185 (non-coding, skips Ex15) | 4.03 | 6.9 × 10⁻⁹ |
| transcript11 (coding, NMD) | 4.94 | 1.6 × 10⁻²³ |
| transcript157 (coding, extra novel intron) | 4.53 | 1.3 × 10⁻⁶ |

All three significant after Bonferroni correction across 3 tests.

## Two views of the same data

`novel_isoforms_per_cluster_UNIQUE_pct.tsv` — **rigorous baseline** (each read assigned to exactly one transcript). Use this for absolute frequencies.

`novel_isoforms_per_cluster_pct.tsv` — re-quant view (a read can match multiple compatible transcripts). Numbers ~10× larger because reads ambiguous between novel and canonical structures count toward both. The neuron-enrichment direction is identical between views.

## Cassette intersection

Only `transcript185` overlaps the lab's 79-bp Ex15 cassette region (mm10 ~chr5:135,869,720-135,869,800):

```
exon 1: 135,867,778-135,868,519
exon 2: 135,868,828-135,869,052
exon 3: 135,869,155-135,869,294    ← upstream of cassette
       [intron: 135,869,295-135,873,077]    ← cassette sits in this intron
exon 4: 135,873,078-135,874,772    ← downstream of cassette
```

`transcript185` is therefore a **novel cassette-skipping non-coding Srrm3 isoform** with a distinct 4-exon 3'-end structure — biologically consistent with the established 0% Ex15 PSI in mature hippocampus, but with different exon boundaries than canonical Srrm3-201.

## Filtered out as artifacts (8 single-exon "novels")

8 additional Srrm3 single-exon novel calls had `perc_A_downstream` ≥ 50% (up to 70%) with downstream sequences like `aaTTACagatagatagatag` — classic intra-priming artifacts (oligo-dT priming on intronic A-stretches). Excluded from the high-confidence set.

## Caveats requiring follow-up validation

1. **Several junctions in the novel transcripts are non-canonical** (GTAC, CGAC, AGGC instead of GT-AG). Could be real (minor U12 spliceosome) or pbmm2 alignment artifacts. **Recommend RT-PCR + Sanger across novel junctions before acting on these calls.**
2. **SQANTI3 v5.4 RTS_junction = `????`** for all junctions — RT-switching detection appears to have not produced a definitive verdict in this version. Manual direct-repeat scanning at junction sequences would close this gap.
3. **No CAGE-peak / polyA-site cross-reference** (would tighten 5' / 3' end confidence). Not blocking but a clean follow-up.
4. **CB recovery rate is 56%** (matching the targeted_psi baseline). Per-cluster counts are conservative by ~2× because the unmatched 44% don't contribute to the per-cell-type breakdown.
5. **transcript185 spans but does not include the Ex15 cassette.** It is a novel cassette-*skipping* isoform, not a novel cassette-*including* isoform.

## Files

| File | What |
|---|---|
| `srrm3_novel_models.gtf` | IGV-loadable GTF: 3 novel + 4 reference Srrm3 transcripts |
| `srrm3_novel_highconf.tsv` | 3-row high-confidence summary table |
| **`novel_isoforms_per_cluster_UNIQUE_counts.tsv`** | **rigorous unique-only counts (baseline)** |
| **`novel_isoforms_per_cluster_UNIQUE_pct.tsv`** | **rigorous unique-only % composition (baseline)** |
| `novel_isoforms_per_cluster_counts.tsv` | re-quant counts (allows multi-counting) |
| `novel_isoforms_per_cluster_pct.tsv` | re-quant % composition |
| `novel_isoforms_per_read.tsv` | per-read CB + cell-type joined table (5.1 MB) |
| `novel_isoforms_per_cluster_per_mouse.tsv` | replicate breakdown by mouse |
| `novel_isoforms_reads_per_100_cells.tsv` | cell-density-normalized counts |
| `novel_isoforms_per_cluster.png` | stacked bar plot (re-quant view) |
| `novel_isoforms_per_cluster_UNIQUE.png` | stacked bar plot (unique-only view) |
| `ISOFORM_COMPOSITION.md` | composition + enrichment markdown |

Adjacent IsoQuant outputs: `../OUT/` (full discovery) and `../../isoquant_targeted_quant/OUT/` (re-quantification with extended annotation).

Adjacent SQANTI3 output: `../../sqanti3/ribostamp_srrm3_classification.txt_tmp`.
EOF

echo "Wrote $SUM"
ls -lh "$RES" | head -25
