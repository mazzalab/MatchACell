<h1 align="center">🍵 MatchACell</h1>

<p align="center">
  <em>Matcha-grade single-cell consensus.</em><br>
  A Snakemake toolkit for single-cell type/state annotation.
</p>

<p align="center">
  <a href="https://github.com/&lt;OWNER&gt;/MatchACell/actions"><img alt="CI" src="https://github.com/&lt;OWNER&gt;/MatchACell/actions/workflows/ci.yml/badge.svg"></a>
  <img alt="Snakemake" src="https://img.shields.io/badge/snakemake-7.32.4-4A7C3C">
  <img alt="Python" src="https://img.shields.io/badge/python-3.10-7AB661">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-7AB661">
</p>

---

> **Status — Step 1 of a multi-step pipeline.**
> MatchACell is being built incrementally. **What ships today is Step 1**: a
> data-driven QC and **Leiden cluster-stability optimizer** that recommends the
> clustering resolution to carry into annotation (the *MatchA Verdict*). The
> longer-term goal is to run **multiple cell annotators in parallel** as
> Snakemake rules and reconcile them into a consensus. Those steps
> (`annotate`, `consensus`) are on the roadmap but **not yet wired**.

## What it does

Given a single-cell or spatial (e.g. **Xenium**) `.h5ad` of post-segmentation
counts, Step 1:

1. **Data-driven QC** — MAD-based outlier detection (no hard-coded thresholds),
   adaptive to the assay (mitochondrial genes and/or Xenium negative-control
   probes).
2. **Preprocess** — normalize → log1p → optional HVG → PCA → neighbours
   (CPU `scanpy`, or GPU `rapids-singlecell`).
3. **Multi-resolution Leiden** clustering across a sweep of resolutions.
4. **Embeddings** (PCA / t-SNE / UMAP) coloured by every resolution.
5. **Bootstrap Jaccard stability** per cluster, per resolution.
6. **The MatchA Verdict** — an automatic recommendation of the finest *stable*
   resolution for annotation, with the clusters to watch.
7. A **self-contained interactive HTML report** with a resolution switcher, plus
   QC plots, a clustree, transition heatmaps, and machine-readable performance
   metrics.

The matcha-green theme applies to the CLI chrome and the HTML report only;
**scientific plots keep their default colour maps**.

## Repository layout

```
MatchACell/
├── run.py                      # matcha-themed Snakemake 7 launcher
├── Snakefile                   # workflow entrypoint
├── config.yaml                 # samples + Step-1 parameters
├── environment.yml             # CPU conda environment (snakemake 7.32.4 pinned)
├── environment-gpu.yml         # GPU (RAPIDS) environment template
├── workflow/
│   ├── rules/cluster_stability.smk
│   ├── envs/matchacell.yaml     # per-rule conda env
│   └── scripts/matchacell_cluster_stability.py   # the Step-1 engine
├── tests/                      # pytest suite (unit + CPU end-to-end)
├── tools/fetch_pbmc3k.py       # grabs the example dataset
├── docs/                       # full documentation
└── .github/workflows/ci.yml    # pytest on CPU
```

## Quickstart

```bash
# 1. Environment (CPU)
conda env create -f environment.yml   # or: mamba env create -f environment.yml
conda activate matchacell

# 2. Example data
python tools/fetch_pbmc3k.py data/pbmc3k.h5ad

# 3. Point config.yaml at an output dir, then dry-run and run
./run.py -w cluster_stability -c config.yaml -q 8 -n   # dry-run
./run.py -w cluster_stability -c config.yaml -q 8      # real run

# See targets
./run.py --list-workflows
```

Outputs land under `<output_dir>/results/<sample>/matchacell/` — start with
`MatchACell_report.html` and `MatchA_Verdict.txt`.

## CPU vs GPU — read this before comparing runs

The CPU and GPU paths are **not drop-in equivalents**. CPU uses
`pynndescent` + `leidenalg`/`igraph`; GPU uses cuML kNN + cuGraph Leiden. They
produce **different cluster counts at the same resolution**, and GPU Jaccard
scores tend to run systematically lower. Treat them as **distinct backends with
distinct result distributions**, not interchangeable acceleration modes. On
small datasets (< ~50k cells) GPU is often *slower* due to many tiny sequential
re-clusterings — the bootstrap auto-routes to CPU below
`--gpu-cell-threshold`. See [`docs/methods.md`](docs/methods.md#cpu-vs-gpu-backend-divergence).

## Documentation

- [`docs/installation.md`](docs/installation.md) — CPU and GPU setup.
- [`docs/pipeline.md`](docs/pipeline.md) — `run.py`, `config.yaml`, the rule, and how to add future annotators.
- [`docs/methods.md`](docs/methods.md) — the science: QC, stability, the Verdict, backend caveats, known limitations.
- [`docs/cli-and-outputs.md`](docs/cli-and-outputs.md) — full CLI reference and an output-by-output map.

## Testing

```bash
pytest -m "not slow"   # fast unit tests
pytest                 # include the CPU end-to-end smoke test
```

## Roadmap

- [x] **Step 1** — QC + Leiden cluster-stability optimizer (the MatchA Verdict)
- [ ] **Step 2** — multiple cell type/state annotators run in parallel as rules
- [ ] **Step 3** — cross-annotator consensus + unified report
- [ ] Reconcile/document the CPU↔GPU backend divergence for end users
- [ ] Seed cuGraph Leiden for run-to-run reproducibility on GPU

## License

MIT — see [`LICENSE`](LICENSE). Update the author fields in `LICENSE` and
`CITATION.cff` before publishing.
