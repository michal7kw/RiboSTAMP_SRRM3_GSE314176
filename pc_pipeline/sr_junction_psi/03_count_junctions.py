#!/usr/bin/env python3
"""Junction-counting PSI calculator (rMATS-style logic).

Counts reads at the three Srrm3 cassette splice junctions in mm10 coordinates.
Works on both bulk and 10x BAMs. For 10x BAMs, optionally aggregates per
cell barcode using the CB tag for downstream per-cluster PSI.

A read is counted as supporting a junction X→Y if its CIGAR has an N op
whose left edge is within ±SR_JUNC_TOL of X AND whose right edge is within
±SR_JUNC_TOL of Y.

PSI is then computed as:
    PSI = (Inc1 + Inc2) / (Inc1 + Inc2 + 2*Skip)

This is the standard rMATS junction-counting formula. Note the factor of 2
on Skip — because the Skip junction is a SINGLE junction event covering
both the cassette's 5' and 3' splice sites, while inclusion has TWO separate
junctions (one each at 5' and 3' of the cassette). Multiplying Skip by 2
keeps inclusion and skip events on the same per-junction basis.

Usage:
    python 03_count_junctions.py --bam INPUT.bam --output OUT.tsv
        [--per-cell]  (10x mode: split counts per CB tag)
        [--cell-list whitelist.tsv]  (only count cells in this list)

Output schema (bulk mode):
    sample, n_inc1, n_inc2, n_skip, psi, psi_lo95, psi_hi95

Output schema (per-cell mode):
    sample, cb, n_inc1, n_inc2, n_skip, psi, psi_lo95, psi_hi95
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pysam
import pandas as pd
from scipy.stats import binomtest

# Pull config from environment (set by 00_config.sh when running)
CHROM        = os.environ.get("SR_CHROM", "chr5")
INC1_LEFT    = int(os.environ.get("SR_INC1_LEFT", 135869294))
INC1_RIGHT   = int(os.environ.get("SR_INC1_RIGHT", 135869720))
INC2_LEFT    = int(os.environ.get("SR_INC2_LEFT", 135869798))
INC2_RIGHT   = int(os.environ.get("SR_INC2_RIGHT", 135873080))
SKIP_LEFT    = INC1_LEFT
SKIP_RIGHT   = INC2_RIGHT
TOL          = int(os.environ.get("SR_JUNC_TOL", 5))
CB_TAG       = os.environ.get("SR_CB_TAG", "CB")

# A query region wide enough to capture any read whose CIGAR contains a
# junction touching the cassette
QUERY_REGION = (CHROM, max(0, SKIP_LEFT - 200), SKIP_RIGHT + 200)


def junction_match(n_left, n_right, target_left, target_right, tol):
    """Does this N op match the target junction within tolerance?"""
    return abs(n_left - target_left) <= tol and abs(n_right - target_right) <= tol


def classify_read(read):
    """Return one of: 'inc1', 'inc2', 'skip', or None.

    Walks the CIGAR ops of the read, looking for N ops that match each
    junction. A single read can support multiple junctions (e.g. a long
    read covering both inc1 and inc2 of an inclusion event); we return
    the first match found, ordered by relevance.

    For short reads (~90 bp), a single read typically supports at most
    ONE junction.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return None
    if read.cigartuples is None:
        return None

    ref_pos = read.reference_start  # 0-based
    matches = []
    for op, length in read.cigartuples:
        # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8
        if op == 3:  # N — splice/skip
            n_left = ref_pos
            n_right = ref_pos + length
            # Check each target junction (mm10, 0-based)
            if junction_match(n_left, n_right, INC1_LEFT, INC1_RIGHT, TOL):
                matches.append("inc1")
            if junction_match(n_left, n_right, INC2_LEFT, INC2_RIGHT, TOL):
                matches.append("inc2")
            if junction_match(n_left, n_right, SKIP_LEFT, SKIP_RIGHT, TOL):
                matches.append("skip")
            ref_pos += length
        elif op in (0, 7, 8, 2):  # M, =, X, D consume reference
            ref_pos += length
        # I and S consume read, not reference

    if not matches:
        return None
    # If both inc and skip match (rare — could happen with chimeric alignments
    # or alignment artifacts), prefer inc since it's a more specific signal.
    # In practice these almost never co-occur in the same read.
    if "inc1" in matches:
        return "inc1"
    if "inc2" in matches:
        return "inc2"
    return "skip"


def wilson_ci(num_successes, num_trials):
    """Wilson 95% binomial CI. Returns (low, high)."""
    if num_trials == 0:
        return (None, None)
    res = binomtest(num_successes, num_trials)
    lo, hi = res.proportion_ci(method="wilson")
    return (lo, hi)


def compute_psi(n_inc1, n_inc2, n_skip):
    """rMATS-style PSI: (Inc1+Inc2) / (Inc1+Inc2 + 2*Skip).

    Returns (psi, ci_lo, ci_hi). psi is None if the denominator is 0.
    """
    inc_total = n_inc1 + n_inc2
    skip_total = n_skip * 2
    n_trials = inc_total + skip_total
    if n_trials == 0:
        return (None, None, None)
    psi = inc_total / n_trials
    lo, hi = wilson_ci(inc_total, n_trials)
    return (psi, lo, hi)


def process_bam_bulk(bam_path: Path, sample_name: str) -> pd.DataFrame:
    """Bulk mode: aggregate junction counts across all reads in the BAM."""
    print(f"  bulk-processing {bam_path}")
    n_inc1 = n_inc2 = n_skip = 0
    n_total = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(*QUERY_REGION):
            n_total += 1
            cls = classify_read(read)
            if cls == "inc1":
                n_inc1 += 1
            elif cls == "inc2":
                n_inc2 += 1
            elif cls == "skip":
                n_skip += 1

    psi, lo, hi = compute_psi(n_inc1, n_inc2, n_skip)
    print(f"    n_total at locus: {n_total:,}, inc1: {n_inc1}, inc2: {n_inc2}, skip: {n_skip}, PSI: {psi}")

    return pd.DataFrame([{
        "sample": sample_name,
        "n_reads_at_locus": n_total,
        "n_inc1": n_inc1,
        "n_inc2": n_inc2,
        "n_skip": n_skip,
        "psi": psi,
        "psi_lo95": lo,
        "psi_hi95": hi,
    }])


def process_bam_per_cell(
    bam_path: Path,
    sample_name: str,
    cell_whitelist: set | None = None,
) -> pd.DataFrame:
    """Per-cell mode: split junction counts by CB tag.

    Reads without a CB tag are skipped. If cell_whitelist is provided, only
    reads whose CB is in the whitelist are counted.
    """
    print(f"  per-cell-processing {bam_path}")
    counts = {}  # cb -> {"inc1": int, "inc2": int, "skip": int}
    n_total = n_with_cb = n_in_whitelist = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(*QUERY_REGION):
            n_total += 1
            try:
                cb = read.get_tag(CB_TAG)
            except KeyError:
                continue
            n_with_cb += 1
            if cell_whitelist is not None and cb not in cell_whitelist:
                continue
            n_in_whitelist += 1

            cls = classify_read(read)
            if cls is None:
                continue
            d = counts.setdefault(cb, {"inc1": 0, "inc2": 0, "skip": 0})
            d[cls] += 1

    print(f"    n_total at locus: {n_total:,}, with CB: {n_with_cb:,}, in whitelist: {n_in_whitelist:,}")
    print(f"    cells with ≥1 informative read: {len(counts):,}")

    rows = []
    for cb, d in counts.items():
        psi, lo, hi = compute_psi(d["inc1"], d["inc2"], d["skip"])
        rows.append({
            "sample": sample_name,
            "cb": cb,
            "n_inc1": d["inc1"],
            "n_inc2": d["inc2"],
            "n_skip": d["skip"],
            "psi": psi,
            "psi_lo95": lo,
            "psi_hi95": hi,
        })
    return pd.DataFrame(rows)


_BC_RE = __import__("re").compile(r"^([ACGT]+)")


def _strip_suffix(s: str) -> str:
    """Strip any non-ACGT suffix (-1, _sampleN, etc) — barcodes are pure ACGT."""
    m = _BC_RE.match(str(s))
    return m.group(1) if m else str(s)


def load_cell_whitelist(path: Path) -> set:
    """Load barcode whitelist (one per line, or first column of TSV).

    Accepts barcodes in any of these forms in the input file:
      AAACCAAAGCAGGTTC          (bare 16-bp — STARsolo BAM CB tag form)
      AAACCAAAGCAGGTTC-1        (10x cellranger -1 suffix)
      AAACCAAAGCAGGTTC_sample1  (long-read pipeline _sampleN suffix)

    For matching against BAM CB tags, we add BOTH the bare form and the
    "-1" form to the set, since different upstream pipelines use either.
    """
    cbs = set()
    with open(path) as f:
        # Skip header if present (first line that doesn't start with ACGT-only)
        first = True
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cb_raw = line.split("\t")[0].split(",")[0]
            bare = _strip_suffix(cb_raw)
            if first:
                first = False
                # If first non-comment line is the literal header word "barcode"
                # (or anything not all-ACGT), skip it.
                if not bare or len(bare) < 8 or not all(c in "ACGT" for c in bare):
                    continue
            if not bare or len(bare) < 8:
                continue
            cbs.add(bare)
            cbs.add(bare + "-1")
    return cbs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path,
                        help="Input BAM (sorted + indexed)")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output TSV path")
    parser.add_argument("--sample-name", default=None,
                        help="Sample name (default: BAM filename stem)")
    parser.add_argument("--per-cell", action="store_true",
                        help="10x mode: split counts per CB tag")
    parser.add_argument("--cell-whitelist", type=Path, default=None,
                        help="Optional whitelist of cell barcodes to include")
    args = parser.parse_args()

    if not args.bam.exists():
        sys.exit(f"BAM not found: {args.bam}")
    sample_name = args.sample_name or args.bam.stem
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print(f"Junction-PSI calculator for sample: {sample_name}")
    print(f"  Region: {CHROM}:{QUERY_REGION[1]:,}-{QUERY_REGION[2]:,}")
    print(f"  Junctions (mm10):")
    print(f"    Inc1: {INC1_LEFT:,} → {INC1_RIGHT:,}")
    print(f"    Inc2: {INC2_LEFT:,} → {INC2_RIGHT:,}")
    print(f"    Skip: {SKIP_LEFT:,} → {SKIP_RIGHT:,}")
    print(f"  Tolerance: ±{TOL} bp")

    if args.per_cell:
        whitelist = None
        if args.cell_whitelist:
            whitelist = load_cell_whitelist(args.cell_whitelist)
            print(f"  Loaded {len(whitelist):,} barcodes from whitelist")
        df = process_bam_per_cell(args.bam, sample_name, whitelist)
    else:
        df = process_bam_bulk(args.bam, sample_name)

    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {len(df):,} rows to {args.output}")


if __name__ == "__main__":
    main()
