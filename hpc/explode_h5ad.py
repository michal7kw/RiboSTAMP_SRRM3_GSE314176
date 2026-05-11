#!/usr/bin/env python3
"""
Explode a Scanpy .h5ad into plain files (obs.tsv, var.tsv, X.mtx, layers)
that the local R scripts can read directly via Matrix::readMM + fread.

Output filename pattern: <output_dir>/<input_basename>__<part>.<ext>
where <input_basename> is the input filename without the .h5ad suffix.
This means multiple h5ad files can be exploded into the same output dir
without clobbering each other — each gets its own namespaced files.

Files emitted per h5ad:
    <basename>__obs.tsv          cell metadata (barcode + obs columns)
    <basename>__var.tsv          feature metadata (feature_id + var columns)
    <basename>__X.mtx            main matrix, transposed to features x cells
    <basename>__layer_<name>.mtx for each layer

Usage:
    python explode_h5ad.py <output_dir> <input.h5ad> [<input2.h5ad> ...]
    python explode_h5ad.py <output_dir> <glob_pattern.h5ad>

Backward-compat: if called with exactly one input (old signature
"explode_h5ad.py <input.h5ad> <output_dir>"), still works.
"""
import sys
import os
import glob
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import io as sio
from scipy import sparse


def explode_one(in_path: Path, out_dir: Path) -> None:
    base = in_path.stem  # filename without .h5ad
    print(f"\n[{base}]  loading...", flush=True)
    a = ad.read_h5ad(in_path)
    print(f"  shape:  {a.shape}  (cells x features)")
    print(f"  obs cols: {list(a.obs.columns)}")
    print(f"  var cols: {list(a.var.columns)[:8]}{'...' if len(a.var.columns) > 8 else ''}")
    print(f"  layers:   {list(a.layers.keys())}")

    # Cell metadata. Insert 'barcode' as the first column for unambiguous
    # joining downstream (the AnnData index can be duplicated).
    obs = a.obs.copy()
    obs.insert(0, "barcode", a.obs_names.astype(str))
    obs.to_csv(out_dir / f"{base}__obs.tsv", sep="\t", index=False)
    print(f"  wrote {base}__obs.tsv  ({len(obs)} rows)")

    # Feature metadata.
    var = a.var.copy()
    var.insert(0, "feature_id", a.var_names.astype(str))
    var.to_csv(out_dir / f"{base}__var.tsv", sep="\t", index=False)
    print(f"  wrote {base}__var.tsv  ({len(var)} rows)")

    # Main X matrix — transpose to features x cells so it matches the 10x
    # MTX convention every R single-cell tool expects (rows=features,
    # cols=cells).
    X = a.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    sio.mmwrite(out_dir / f"{base}__X.mtx", X.T)
    print(f"  wrote {base}__X.mtx  (features x cells)")

    # All layers, same orientation.
    for layer_name in a.layers.keys():
        L = a.layers[layer_name]
        if not sparse.issparse(L):
            L = sparse.csr_matrix(L)
        out_file = out_dir / f"{base}__layer_{layer_name}.mtx"
        sio.mmwrite(out_file, L.T)
        print(f"  wrote {out_file.name}")


def main() -> None:
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)

    # Backward-compat: old signature was (input.h5ad, output_dir). If first
    # arg ends with .h5ad and second is a directory (or doesn't exist as a
    # file), assume the old order.
    if (args[0].endswith(".h5ad") and not args[1].endswith(".h5ad")):
        in_paths = [Path(args[0])]
        out_dir = Path(args[1])
    else:
        out_dir = Path(args[0])
        # Expand globs in remaining args
        in_paths = []
        for p in args[1:]:
            if any(c in p for c in "*?["):
                matched = sorted(Path(x) for x in glob.glob(p))
                if not matched:
                    sys.exit(f"glob matched nothing: {p}")
                in_paths.extend(matched)
            else:
                in_paths.append(Path(p))

    out_dir.mkdir(parents=True, exist_ok=True)
    for p in in_paths:
        if not p.is_file():
            sys.exit(f"input not found: {p}")
        explode_one(p, out_dir)


if __name__ == "__main__":
    main()
