#!/usr/bin/env python3
"""Build a per-mouse STARsolo CB whitelist from the authors' obs.tsv.

We use the authors' filtered cell barcodes as the STARsolo whitelist
(rather than the 10x master 3M-february-2018 whitelist) because:

  1. We only care about cells the authors retained — those are the only
     ones with cluster-type labels we can use downstream.
  2. STARsolo will Hamming-1 error-correct against the supplied whitelist,
     so passing only ~6,500 real cells gives high-precision CB recovery.
  3. We don't have to download or license the 10x master whitelist.

Output: one line per 16-bp barcode (no "-1" suffix, no header).
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

BC_RE = re.compile(r"^([ACGT]+)")


def strip_suffix(s: str) -> str:
    """Strip any non-ACGT suffix (-1, _sampleN, etc) — barcodes are pure ACGT."""
    m = BC_RE.match(str(s))
    return m.group(1) if m else str(s)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--obs", required=True, type=Path,
                   help="Short-read obs.tsv (must have 'barcode' column)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output whitelist .txt (one barcode per line)")
    args = p.parse_args()

    if not args.obs.exists():
        sys.exit(f"ERROR: obs.tsv not found: {args.obs}")

    df = pd.read_csv(args.obs, sep="\t")
    if "barcode" not in df.columns:
        sys.exit(f"ERROR: 'barcode' column missing in {args.obs}")

    bare = df["barcode"].astype(str).map(strip_suffix)
    # Drop any non-16bp entries (defensive — should not occur for 10x v3)
    keep_mask = bare.str.len() == 16
    if (~keep_mask).any():
        n_drop = (~keep_mask).sum()
        print(f"  WARN: dropping {n_drop} barcodes whose stripped length != 16",
              file=sys.stderr)
    bare = bare[keep_mask].drop_duplicates()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    bare.to_csv(args.output, header=False, index=False)
    print(f"Wrote {len(bare):,} unique 16-bp barcodes → {args.output}")


if __name__ == "__main__":
    main()
