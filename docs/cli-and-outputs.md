# CLI reference & outputs

## Running the engine directly

The Snakemake rule calls the engine, but you can also run it standalone:

```bash
python workflow/scripts/matchacell_cluster_stability.py \
    --input data/pbmc3k.h5ad \
    --outdir results/pbmc3k/matchacell \
    --backend auto \
    --n-iter 1000
```

## All flags

### Required
| Flag | Default | Description |
| --- | --- | --- |
| `--input` | — | Input `.h5ad` (raw-ish counts). |
| `--outdir` | — | Output directory (created if absent). |

### Backend & performance
| Flag | Default | Description |
| --- | --- | --- |
| `--backend` | `auto` | `auto` \| `gpu` \| `cpu` for the main pipeline. |
| `--stability-backend` | `auto` | Backend for the bootstrap loop specifically. |
| `--gpu-cell-threshold` | `50000` | With `auto`, datasets smaller than this run the bootstrap on CPU. |
| `--rmm-pool-fraction` | `0.5` | Fraction of free VRAM pre-allocated for the RMM pool. |
| `--rmm-managed` | off | Use managed (unified) memory for the RMM pool. |
| `--skip-perf` | off | Disable performance profiling outputs. |
| `--no-color` | off | Disable matcha CLI colours. |

### Clustering
| Flag | Default | Description |
| --- | --- | --- |
| `--resolutions` | `0.1 0.2 0.3 0.5 0.8 1.0 1.5 2.0` | Leiden resolution sweep. |
| `--n-pcs` | `50` | Principal components for PCA/neighbours. |
| `--n-neighbors` | `15` | kNN graph neighbours. |
| `--use-hvg` | `auto` | `auto` \| `yes` \| `no` (auto = HVG when > 2,000 genes). |

### Embeddings
| Flag | Default | Description |
| --- | --- | --- |
| `--skip-embeddings` | off | Skip embedding plots. |
| `--skip-tsne` | off | Skip t-SNE (the slowest embedding). |
| `--skip-umap` | off | Skip UMAP. |

### QC
| Flag | Default | Description |
| --- | --- | --- |
| `--qc-nmads` | `5.0` | MAD sensitivity for count/gene metrics. |
| `--qc-nmads-pct` | `3.0` | MAD sensitivity for percentage metrics. |
| `--min-counts` | `10` | Absolute floor on total counts. |
| `--min-genes` | `5` | Absolute floor on genes per cell. |
| `--mt-pct-hard` | `20.0` | Hard mitochondrial-percent ceiling. |
| `--skip-qc` | off | Skip QC entirely. |

### Stability, transitions, report
| Flag | Default | Description |
| --- | --- | --- |
| `--n-iter` | `100` | Bootstrap iterations per resolution. |
| `--fraction` | `0.8` | Subsample fraction per iteration. |
| `--skip-stability` | off | Skip the bootstrap (and the Verdict). |
| `--skip-transitions` | off | Skip transition heatmaps. |
| `--skip-html` | off | Skip the interactive HTML report. |

### Verdict thresholds
| Flag | Default | Description |
| --- | --- | --- |
| `--verdict-high` | `0.85` | "High stability" reference line. |
| `--verdict-risk` | `0.60` | Risk threshold for unstable clusters. |
| `--verdict-min-median` | `0.75` | Minimum median stability to be acceptable. |
| `--verdict-max-risk-frac` | `0.25` | Max fraction of risky clusters allowed. |

### Misc
| Flag | Default | Description |
| --- | --- | --- |
| `--edge-weight-threshold` | `0.05` | Clustree edge-weight cutoff. |
| `--seed` | `0` | RNG seed (CPU bootstrap is reproducible; see GPU caveat in methods). |

> Note: when set through the Snakemake rule, the config's `n_iter` defaults to
> `1000`; the script's own default (standalone) is `100`.

## Output map

All paths are relative to `--outdir`
(`<output_dir>/results/<sample>/matchacell/` under the workflow).

```
clustered_multi_resolution.h5ad     # tracked output: all leiden_<res> + embeddings + uns
MatchA_Verdict.txt                  # the recommendation, human-readable
MatchACell_report.html              # self-contained interactive report (start here)
clustree_<thr>.png                  # cluster tree across resolutions (if pyclustree)

qc/
  qc_summary.csv
  qc_violin_kept_vs_removed.png
  qc_hist_mad_bounds.png
  qc_counts_vs_genes.png
  qc_barcode_rank.png
  qc_removal_by_criterion.png
  qc_highest_expr_genes.png

embeddings/
  embeddings_leiden_<res>.png       # PCA/t-SNE/UMAP coloured by each resolution
  embeddings_qc_metrics.png

stability/
  Cluster_Stability_Summary.xlsx    # one sheet per resolution
  global_stability_comparison_boxplot.png
  resolution_summary.csv

transitions/
  <parent>_to_<child>_heatmap_counts.png
  All_cluster_transitions.xlsx

performance/
  performance_metrics.csv           # per-stage wall-clock
  stability_timing.csv              # per-resolution throughput
  performance_summary.json          # run metadata (backend, GPU, cells, total time)
  stage_durations.png
  stability_time_per_resolution.png
```

Start with `MatchACell_report.html` (it bundles QC, stability, an embedding with
a resolution switcher, transitions, and performance) and `MatchA_Verdict.txt`.
