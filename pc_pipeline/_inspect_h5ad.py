#!/usr/bin/env python3
"""One-off inspection of GSE314176 h5ad files. Prints layers, X dtype,
shape, obs+var columns, and counts of unique cell-type labels."""
import sys
import glob
import anndata as ad
import pandas as pd

paths = sys.argv[1:] or sorted(glob.glob(
    "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/processed/GSE314176_RAW/*.h5ad"))
if not paths:
    print("No h5ad files found", file=sys.stderr); sys.exit(1)

for p in paths:
    print(f"\n{'='*70}\n  {p}\n{'='*70}")
    a = ad.read_h5ad(p, backed="r")
    print(f"  shape (cells × features): {a.shape}")
    print(f"  X.dtype: {getattr(a.X, 'dtype', type(a.X))}")
    print(f"  layers (name -> dtype):")
    for name in a.layers.keys():
        L = a.layers[name]
        print(f"    {name}: {getattr(L, 'dtype', type(L))}")
    print(f"  obs columns ({len(a.obs.columns)}):")
    for c in a.obs.columns:
        print(f"    {c!r}")
    print(f"  var columns ({len(a.var.columns)}):")
    for c in a.var.columns[:10]:
        print(f"    {c!r}")
    if len(a.var.columns) > 10:
        print(f"    ... and {len(a.var.columns)-10} more")
    # Cell-type column candidates
    for col in a.obs.columns:
        if any(k in col.lower() for k in ("type", "assign", "cluster",
                                          "annotation", "label", "celltype")):
            print(f"\n  unique values in {col!r}:")
            vc = a.obs[col].value_counts(dropna=False)
            for label, n in vc.items():
                print(f"    {n:>5}  {label!r}")
