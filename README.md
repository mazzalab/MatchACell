<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/logo/matchacell-logo-dark.svg">
    <img src="docs/logo/matchacell-logo.svg" alt="MatchACell" width="340">
  </picture>
</p>

<p align="center">
  A Snakemake toolkit for single-cell type/state annotation.
</p>

<p align="center">
  <a href="https://github.com/mazzalab/MatchACell/actions"><img alt="CI" src="https://github.com/mazzalab/MatchACell/actions/workflows/ci.yml/badge.svg"></a>
  <img alt="Snakemake" src="https://img.shields.io/badge/snakemake-7.32.4-4A7C3C">
  <img alt="Python" src="https://img.shields.io/badge/python-3.10-7AB661">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-7AB661">
</p>

---

> **Status: Steps 1 and 2 are available.** Step 1 cleans and clusters the data
> and recommends a clustering resolution, the *MatchA Verdict*. Step 2 annotates
> that clustering with up to eight independent annotators. Cross-annotator
> consensus (Step 3) is on the roadmap.

## How it works

```
input .h5ad (or Seurat .rds)
  │
  ├─ Step 0  stage every sample as .h5ad
  ├─ Step 1  QC → multi-resolution Leiden → bootstrap stability → MatchA Verdict
  ├─ Step 2  annotators, in parallel, at the Verdict resolution
  └─ Step 3  consensus across annotators (planned)
```

### Step 1: clustering

Given a single-cell or spatial (e.g. **Xenium**) `.h5ad` of post-segmentation
counts, Step 1 runs:

1. **Data-driven QC**: MAD-based outlier detection (no hard-coded thresholds),
   adaptive to the assay (mitochondrial genes and/or Xenium negative-control
   probes).
2. **Preprocessing**: normalize → log1p → optional HVG → PCA → neighbours
   (CPU `scanpy`, or GPU `rapids-singlecell`).
3. **Multi-resolution Leiden** clustering across a sweep of resolutions.
4. **Embeddings** (PCA / t-SNE / UMAP) coloured by every resolution.
5. **Bootstrap Jaccard stability** per cluster, per resolution.
6. **The MatchA Verdict**: the finest *stable* resolution for annotation, with
   the clusters to watch.
7. A **self-contained interactive HTML report** with a resolution switcher, plus
   QC plots, a clustree, transition heatmaps, and performance metrics.

### Step 2: annotation

The `annotation` target runs every enabled annotator as its own Snakemake rule,
with its own conda environment, on the clustering the Verdict recommends.

| Annotator | Approach | Labels | Runs when |
| --- | --- | --- | --- |
| [ScoreGenes](docs/annotators/score-genes.md) | Marker-gene signatures, Scanpy `score_genes` | per cluster | always |
| [AddModuleScore](docs/annotators/addmodulescore.md) | Marker-gene signatures, Seurat `AddModuleScore` (R) | per cluster | always |
| [CIA](docs/annotators/cia.md) | Marker-gene signatures, cells classified by CIA | per cell, summarised per cluster | always |
| [CellTypist](docs/annotators/celltypist.md) | Pretrained CellTypist models | per cell, majority vote per cluster | always |
| [scANVI](docs/annotators/scanvi.md) | Reference mapping with scvi-tools | per cell | `scanvi.reference_file` is set |
| [scParadise](docs/annotators/scparadise.md) | Pretrained scAdam model | per cell | `scparadise.model_dir` is set |
| [CyteType](docs/annotators/cytetype.md) | CyteType AI annotation service | per cluster | a CyteType API token is set |
| [CellTypeAI](docs/annotators/celltypeai.md) | Local LLM ensemble via Ollama | per cluster | `celltypeai.tissue` is set |

The three signature-based annotators share one Excel workbook of marker genes,
`matchacell_annotation.annot_file`. The [annotators overview](docs/annotators/README.md)
explains how to prepare it, what each annotator needs, and where its results go.

## Quickstart

```bash
# 1. Launcher environment (CPU)
conda env create -f environment.yml   # or: mamba env create -f environment.yml
conda activate matchacell

# 2. Example data
python tools/fetch_pbmc3k.py data/pbmc3k.h5ad

# 3. In config.yaml, set output_dir, samples and, for Step 2,
#    matchacell_annotation.annot_file (your marker-gene workbook).

# 4. Step 1: QC, clustering and the MatchA Verdict
./run.py -w cluster_stability -c config.yaml -q 8 -n   # dry-run
./run.py -w cluster_stability -c config.yaml -q 8

# 5. Step 2: every enabled annotator (runs Step 1 first if needed)
./run.py -w annotation -c config.yaml -q 8

# On a PBS cluster, submit every job with qsub instead (see docs/pipeline.md)
./run.py -w annotation -c config.yaml -q 1 -cl -qu workq -j 20

# List the targets
./run.py --list-workflows
```

Each rule's conda environment is built on first use, so the first run of each
step takes longer. See [`docs/installation.md`](docs/installation.md).

> **Before running Step 2, check two defaults in `config.yaml`:**
>
> - Step 1 keeps only the 2,000 most variable genes in the object the
>   annotators read. Set `matchacell_cluster_stability.extra: "--use-hvg no"` so
>   pretrained models and signatures keep their genes
>   ([why](docs/annotators/README.md#keep-every-gene-for-annotation)).
> - CellTypeAI is on by default (`celltypeai.tissue: "Immune system"`) and needs
>   a running Ollama server. Set `tissue: ""` if you don't have one.

## Outputs

```
<output_dir>/results/<sample>/
├── h5/<sample>.h5ad                        # staged input
└── matchacell/
    ├── MatchACell_report.html              # Step 1 report: start here
    ├── MatchA_Verdict.txt                  # the recommended resolution
    ├── clustered_multi_resolution.h5ad     # every leiden_<res> column + embeddings
    ├── clustered_multi_resolution.rds      # Seurat companion (read by AddModuleScore)
    ├── qc/  embeddings/  stability/  transitions/  performance/
    └── annotation/
        ├── ScoreGenes/  AddModuleScore/  CIA/  CellTypist/
        └── scANVI/  scParadise/  CyteType/  CellTypeAI/     # when enabled
```

[`docs/cli-and-outputs.md`](docs/cli-and-outputs.md) maps every Step 1 file;
each annotator's page lists its own outputs.

## Repository layout

```
MatchACell/
├── run.py                    # matcha-themed Snakemake 7 launcher
├── Snakefile                 # entrypoint: includes every rule, defines the targets
├── config.yaml               # samples, Step 1 parameters, one section per annotator
├── environment.yml           # launcher environment (CPU, snakemake 7.32.4 pinned)
├── environment-gpu.yml       # GPU (RAPIDS) environment template
├── workflow/
│   ├── rules/                # one rule module per step and per annotator
│   ├── scripts/              # what the rules run (Python and R)
│   └── envs/                 # one conda environment per rule
├── tests/                    # pytest suite (unit, scheduling, CPU end-to-end)
├── tools/fetch_pbmc3k.py     # grabs the example dataset
├── docs/                     # documentation; docs/annotators/ has one page per annotator
└── .github/workflows/ci.yml  # pytest on CPU
```

The matcha-green theme applies to the CLI chrome and the HTML report only;
**scientific plots keep their default colour maps**.

## CPU vs GPU: read this before comparing runs

The CPU and GPU paths are **not drop-in equivalents**. CPU uses
`pynndescent` + `leidenalg`/`igraph`; GPU uses cuML kNN + cuGraph Leiden. They
produce **different cluster counts at the same resolution**, and GPU Jaccard
scores tend to run systematically lower. Treat them as **distinct backends with
distinct result distributions**, not interchangeable acceleration modes. On
small datasets (< ~50k cells) GPU is often *slower* due to many tiny sequential
re-clusterings; the bootstrap auto-routes to CPU below
`--gpu-cell-threshold`. See [`docs/methods.md`](docs/methods.md#cpu-vs-gpu-backend-divergence).

## Documentation

- [`docs/installation.md`](docs/installation.md): the launcher environment, GPU
  setup, and what each annotator needs to run.
- [`docs/pipeline.md`](docs/pipeline.md): `run.py`, `config.yaml`, targets and
  rules, and how to add an annotator.
- [`docs/methods.md`](docs/methods.md): the Step 1 science: QC, stability, the
  Verdict, backend caveats, known limitations.
- [`docs/cli-and-outputs.md`](docs/cli-and-outputs.md): Step 1 engine flags and
  an output-by-output map.
- [`docs/annotators/`](docs/annotators/README.md): Step 2 conventions, the
  signature workbook, and one page per annotator.

## Testing

```bash
pytest -m "not slow"   # fast unit tests
pytest                 # include the CPU end-to-end smoke test
```

## Roadmap

- [x] **Step 1**: QC + Leiden cluster-stability optimizer (the MatchA Verdict)
- [x] **Step 2**: eight annotators run in parallel as rules (ScoreGenes,
  AddModuleScore, CIA, CellTypist, scANVI, scParadise, CyteType, CellTypeAI)
- [ ] **Step 3**: cross-annotator consensus + unified report
- [ ] Reconcile/document the CPU↔GPU backend divergence for end users
- [ ] Seed cuGraph Leiden for run-to-run reproducibility on GPU

## License

MIT; see [`LICENSE`](LICENSE). Update the author fields in `LICENSE` and
`CITATION.cff` before publishing.
