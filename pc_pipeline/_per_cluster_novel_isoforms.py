#!/usr/bin/env python3
"""
Per-cell-type expression of novel Srrm3 isoforms.

Joins:
  IsoQuant read_assignments  (read_id -> isoform_id)
  targeted_psi per_read.tsv  (read_id -> cb_matched, sample/srr)
  authors' obs.tsv           (sample, barcode -> Cell Assignments Grouped)

For each novel Srrm3 transcript, output per-cluster + per-mouse counts.

Run with the same conda env as targeted_psi/04 (anything with pandas).
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths (mirror targeted_psi/04 conventions)
# ---------------------------------------------------------------------------

BASE_DIR = os.environ.get("BASE_DIR", "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data")
ISOQ_OUT = Path(BASE_DIR) / "isoquant_targeted_quant" / "OUT"
TGT_PERREAD = Path(BASE_DIR) / "targeted_psi" / "per_read" / "all_samples.per_read.tsv"
META_DIR = Path(BASE_DIR) / "metadata"
RES_DIR = Path(BASE_DIR) / "isoquant_targeted" / "results"
RES_DIR.mkdir(parents=True, exist_ok=True)

# SRR ↔ author mouse GSM (from targeted_psi/04)
SRR_TO_GSM = {
    "SRR36480452": "GSM9380801",  # mouse3
    "SRR36480453": "GSM9380800",  # mouse2
    "SRR36480454": "GSM9380799",  # mouse1
}

# The 3 high-confidence novel + 1 canonical reference for context
TRANSCRIPTS_OF_INTEREST = {
    "transcript11.chr5.nnic":  "novel — 11 exons, terminal exon shift, alt polyA",
    "transcript157.chr5.nnic": "novel — extra novel intron, exon elongation",
    "transcript185.chr5.nnic": "novel — 4-exon 3'-end form, skips Ex15 region",
    "ENSMUST00000005077.6":    "canonical Srrm3-201 (most-expressed reference)",
    "ENSMUST00000144211.1":    "Srrm3-204 reference (parent of transcript185)",
}


def load_read_assignments() -> pd.DataFrame:
    """Read IsoQuant read_assignments.tsv.gz, keep read_id + isoform_id."""
    p = ISOQ_OUT / "OUT.read_assignments.tsv.gz"
    rows = []
    with gzip.open(p, "rt") as fh:
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
            rows.append((parts[ri], parts[ii], parts[ai]))
    df = pd.DataFrame(rows, columns=["read_id", "isoform_id", "assignment_type"])
    return df


def load_per_read_cb() -> pd.DataFrame:
    """Read targeted_psi all_samples.per_read.tsv, keep read_id + cb_matched + sample."""
    df = pd.read_csv(TGT_PERREAD, sep="\t",
                     dtype={"cb_matched": "string", "sample": "string"})
    df["cb_matched"] = df["cb_matched"].fillna("")
    return df[["sample", "read_id", "cb_matched", "cb_orientation"]]


def load_author_obs() -> pd.DataFrame:
    parts = []
    for srr, gsm in SRR_TO_GSM.items():
        path = next(META_DIR.glob(f"{gsm}_longread_normed_counts_*__obs.tsv"), None)
        if path is None:
            raise FileNotFoundError(f"obs.tsv missing for {gsm}")
        df = pd.read_csv(path, sep="\t")
        df["mouse_gsm"] = gsm
        df["srr"] = srr
        df["bare_barcode"] = df["barcode"].astype(str).str.split("_", n=1).str[0]
        parts.append(df[["srr", "mouse_gsm", "bare_barcode",
                          "Cell Assignment", "Cell Assignments Grouped"]])
    out = pd.concat(parts, ignore_index=True)
    return out.rename(columns={
        "Cell Assignment": "cell_assignment",
        "Cell Assignments Grouped": "cell_type",
    })


def main() -> None:
    print("Loading IsoQuant read_assignments...")
    ra = load_read_assignments()
    print(f"  {len(ra):,} read assignments total")

    # Filter to transcripts of interest
    ra_subset = ra[ra["isoform_id"].isin(TRANSCRIPTS_OF_INTEREST.keys())].copy()
    print(f"  {len(ra_subset):,} reads on transcripts of interest")
    print("    breakdown:")
    print(ra_subset.groupby("isoform_id").size().to_string())

    print("\nLoading targeted_psi per_read CBs...")
    pr = load_per_read_cb()
    print(f"  {len(pr):,} per-read rows")

    # Join read_assignments × per_read on read_id
    joined = ra_subset.merge(pr, on="read_id", how="left", indicator="cb_match")
    cb_recovered = (joined["cb_matched"] != "").sum()
    print(f"  {cb_recovered:,} of {len(joined):,} novel-iso reads have a CB ({100*cb_recovered/max(1,len(joined)):.1f}%)")

    # Diagnostics on missing CB
    missing = joined[joined["cb_matched"] == ""]
    print(f"  reads with no CB: {len(missing):,}")
    if len(missing):
        print("    by cb_match indicator:")
        print(missing.groupby("cb_match", observed=False).size().to_string())

    print("\nLoading author obs (long-read library cell types)...")
    obs = load_author_obs()
    print(f"  {len(obs):,} cells across 3 mice")

    # Join (srr, cb) -> cell type
    joined = joined.rename(columns={"sample": "srr"})
    final = joined.merge(
        obs, left_on=["srr", "cb_matched"],
        right_on=["srr", "bare_barcode"], how="left",
    )
    final["cell_type"] = final["cell_type"].fillna("UNASSIGNED")

    # Per-transcript per-cell-type counts
    print("\n=== Per-cell-type read counts ===")
    pivot = (final.groupby(["isoform_id", "cell_type"]).size()
                  .unstack(fill_value=0))
    pivot["TOTAL"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("TOTAL", ascending=False)
    print(pivot.to_string())

    # Save full per-cluster table
    out_full = RES_DIR / "novel_isoforms_per_cluster.tsv"
    pivot.to_csv(out_full, sep="\t")
    print(f"\nWrote {out_full}")

    # Also: per (transcript, mouse, cell_type) for replicate detail
    print("\n=== Per-mouse × per-cell-type detail ===")
    detail = (final.groupby(["isoform_id", "srr", "cell_type"]).size()
                   .unstack(fill_value=0))
    detail["TOTAL"] = detail.sum(axis=1)
    print(detail.to_string())

    detail_path = RES_DIR / "novel_isoforms_per_cluster_per_mouse.tsv"
    detail.to_csv(detail_path, sep="\t")
    print(f"\nWrote {detail_path}")

    # Also save just the joined per-read table for downstream slicing
    keep = final[["read_id", "isoform_id", "assignment_type",
                  "srr", "cb_matched", "cb_orientation",
                  "cell_assignment", "cell_type", "mouse_gsm"]]
    keep_path = RES_DIR / "novel_isoforms_per_read.tsv"
    keep.to_csv(keep_path, sep="\t", index=False)
    print(f"Wrote {keep_path}")


if __name__ == "__main__":
    main()
