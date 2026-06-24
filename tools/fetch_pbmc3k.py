#!/usr/bin/env python
"""Download the pbmc3k example dataset and save it as data/pbmc3k.h5ad.

Usage:
    python tools/fetch_pbmc3k.py [output_path]
"""
import sys
from pathlib import Path

import scanpy as sc


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/pbmc3k.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata = sc.datasets.pbmc3k()  # ~2,700 cells, raw counts
    adata.write(out)
    print(f"Wrote {adata.n_obs} cells x {adata.n_vars} genes -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
