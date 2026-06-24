# Installation

MatchACell is a Snakemake workflow. The launcher (`run.py`) drives the
**Snakemake 7.x** Python API and is pinned to `7.32.4`; the 8.x API is
incompatible.

## CPU (recommended starting point)

```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate matchacell
```

This installs Python 3.10, `snakemake-minimal=7.32.4`, the single-cell stack
(`scanpy`, `anndata`, `leidenalg`, `python-igraph`, `pynndescent`,
`umap-learn`, `scikit-learn`), plotting/report dependencies (`matplotlib`,
`seaborn`, `plotly`, `openpyxl`), and `pytest`. `pyclustree` is installed via
pip and is optional — the pipeline degrades gracefully if it is absent.

Verify:

```bash
python -c "import snakemake, scanpy; print(snakemake.__version__, scanpy.__version__)"
./run.py --list-workflows
```

## GPU (optional, RAPIDS)

GPU acceleration uses [`rapids-singlecell`](https://rapids-singlecell.readthedocs.io)
on top of RAPIDS (`cuml`, `cugraph`, `rmm`, `cupy`). These libraries are
**CUDA-version- and platform-specific**, so `environment-gpu.yml` is a
**template** — edit `cuda-version` and the RAPIDS pins to a single coherent
RAPIDS release that matches your driver/toolkit before creating the env.

```bash
# Edit environment-gpu.yml first (cuda-version, cuml/cugraph/rmm/cupy)
conda env create -f environment-gpu.yml
conda activate matchacell-gpu
```

Then either set `backend: "auto"` (uses GPU when available) or `backend: "gpu"`
in `config.yaml`, or pass `--backend gpu` via the rule's `extra` flags.

> **Before you benchmark GPU vs CPU**, read the backend-divergence section in
> [`methods.md`](methods.md#cpu-vs-gpu-backend-divergence). The two backends are
> not numerically equivalent, and GPU is frequently *slower* on small datasets.

## Per-rule conda environments

`workflow/rules/cluster_stability.smk` declares
`conda: "../envs/matchacell.yaml"`. When you run with conda integration enabled
(the default; disable with `--no-conda`), Snakemake builds that environment per
rule. If you would rather reuse a prebuilt environment, point the `conda:`
directive at its path, or run with `--no-conda` inside an already-activated
`matchacell` env.

## Example dataset

```bash
python tools/fetch_pbmc3k.py data/pbmc3k.h5ad
```

This downloads the ~2,700-cell pbmc3k benchmark and writes it to
`data/pbmc3k.h5ad`, matching the default `samples:` entry in `config.yaml`.
