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

Each rule runs in its own conda environment, which Snakemake builds on first use
when conda integration is enabled (the default; disable with `--no-conda`). The
first run of each step therefore spends a while creating environments.

| Environment | Rules | Contents |
| --- | --- | --- |
| `matchacell.yaml` | `matchacell_cluster_stability`, `score_genes` | Python 3.10 single-cell stack |
| `addmodulescore.yaml` | `stage_h5`, `matchacell_cluster_stability_rds`, `addmodulescore` | R 4.3, Seurat ≥ 5, zellkonverter, rhdf5 |
| `cia.yaml` | `cia` | Python 3.10 + `cia-python` |
| `celltypist.yaml` | `celltypist` | Python 3.10 + `celltypist` |
| `scanvi.yaml` | `scanvi` | Python 3.10 + `scvi-tools`, PyTorch, JAX |
| `scparadise.yaml` | `scparadise` | Python 3.10 + `scparadise==1.1.0`, PyTorch |
| `cytetype.yaml` | `cytetype` | Python 3.12 + `cytetype` |
| `celltypeai.yaml` | `celltypeai` | Python 3.12 + `celltypeai`, `ollama` |

To reuse a prebuilt environment instead, point that rule's `conda:` directive at
its path, or run with `--no-conda` inside an environment that has everything the
rules need.

### What annotators need besides their environment

| Annotator | Also needs |
| --- | --- |
| ScoreGenes, AddModuleScore, CIA | a [signature workbook](annotators/README.md#the-signature-workbook) (`matchacell_annotation.annot_file`) |
| CellTypist | internet access on the first run, to download its models |
| scANVI | an annotated reference `.h5ad` |
| scParadise | a downloaded scAdam model ([setup](annotators/scparadise.md#model-and-environment)) |
| CyteType | a CyteType API token and internet access |
| CellTypeAI | a running Ollama server with the model pulled ([setup](annotators/celltypeai.md#setting-up-ollama)) |

### mamba 2 and Snakemake 7

Snakemake 7 creates environments with `mamba` by default. With mamba 2.x
(observed with 2.8.1), environment creation fails with:

```
CreateCondaEnvironmentException: ...
error    libmamba Non-conda folder exists at prefix - aborting.
```

`run.py` avoids this by creating environments with conda by default
(`--conda-frontend conda`); pass `--conda-frontend mamba` only with mamba 1.x.
When calling Snakemake directly, add `--conda-frontend conda` yourself.

## Example dataset

```bash
python tools/fetch_pbmc3k.py data/pbmc3k.h5ad
```

This downloads the ~2,700-cell pbmc3k benchmark and writes it to
`data/pbmc3k.h5ad`, matching the default `samples:` entry in `config.yaml`.
