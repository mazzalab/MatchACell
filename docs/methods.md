# Methods & scientific notes

This document describes what the Step-1 engine
(`workflow/scripts/matchacell_cluster_stability.py`) computes and how to read
its output, including the important caveats around the two compute backends.

## 1. Data-driven QC

Cells are filtered using **median absolute deviation (MAD)** outlier detection
rather than fixed thresholds, so the criteria adapt to each dataset.

For a metric *x* with median *m* and `MAD`, a cell is flagged when

```
|x - m| > nmads · MAD
```

(both tails). When `MAD = 0` no cell is flagged. Two sensitivities are used:
`--qc-nmads` (default 5) for count- and gene-based metrics, and the stricter
`--qc-nmads-pct` (default 3) for percentage metrics.

Metrics considered:

- `log1p_total_counts`, `log1p_n_genes_by_counts`, and the top-genes fraction
  (`pct_counts_in_top_20`) — flagged with `--qc-nmads`.
- `pct_counts_mt` (if mitochondrial genes are detected) — flagged with the
  stricter percentage sensitivity *and* a hard ceiling `--mt-pct-hard`
  (default 20%).
- `pct_counts_control` (if Xenium negative-control probes are detected) —
  flagged with the stricter percentage sensitivity.
- A small absolute floor (`--min-counts`, `--min-genes`) removes obviously empty
  barcodes.

**Feature classes.** Mitochondrial genes are matched by the `MT-`/`MT.` prefix.
Control/probe features are matched against Xenium-style patterns
(`negcontrolprobe`, `negcontrolcodeword`, `blank`, `antisense`,
`unassignedcodeword`, `deprecatedcodeword`, `intergenic`, `genomic`). Detected
control features are **dropped from the expression matrix** after QC so they do
not influence PCA or clustering.

QC diagnostics (violins, MAD-bound histograms, counts-vs-genes scatters, a
barcode-rank knee plot, a per-criterion removal bar chart, and the
highest-expressed-genes plot) are written to `qc/`.

## 2. Preprocessing, clustering, embeddings

`normalize_total → log1p → (optional HVG, top 2,000 genes) → PCA → neighbours`.
HVG selection is automatic when the matrix has more than 2,000 genes
(`--use-hvg auto`). Leiden is then run across the resolution sweep
(`--resolutions`, default `0.1 … 2.0`), producing one `leiden_<res>` column per
resolution. PCA / t-SNE / UMAP embeddings are coloured by every resolution and
saved to `embeddings/`. A `clustree` is written when `pyclustree` is installed.

## 3. Bootstrap Jaccard stability

For each resolution the engine estimates how reproducible each cluster is under
resampling:

1. Subsample a fraction (`--fraction`, default 0.8) of cells without
   replacement.
2. Re-run neighbours + Leiden on the subsample's PCA coordinates.
3. Build a contingency table between original and re-clustered labels and take,
   for each original cluster, the **maximum Jaccard index** against any new
   cluster.
4. Average that best-match Jaccard over `--n-iter` iterations.

The result is a per-cluster stability score in `[0, 1]`. A boxplot across
resolutions and a per-resolution Excel summary are written to `stability/`.

On GPU the full PCA matrix is staged to device memory once and every iteration
gathers its subsample on-device, avoiding thousands of host→device copies.

## 4. The MatchA Verdict

The Verdict recommends a resolution to carry into annotation.

A resolution is **acceptable** when its **median** stability ≥ `--verdict-min-median`
(default 0.75) **and** its fraction of risky clusters (< `--verdict-risk`,
default 0.60) ≤ `--verdict-max-risk-frac` (default 0.25). Among acceptable
resolutions, the Verdict picks the **finest** (most clusters) — finer
granularity gives more cell-type/state resolution, provided it stays stable. If
none qualify, it falls back to the resolution with the **highest median
stability**.

The Verdict also lists clusters below the risk threshold to watch (candidates
for merging, marker re-examination, or sub-clustering). It is written to
`MatchA_Verdict.txt`, surfaced in the HTML report, and stored in
`adata.uns["MatchACell"]["verdict_resolution"]`.

Reference thresholds used in plots and the report: **≥ 0.85** "high stability",
**< 0.60** "risk".

## 5. Transitions

Between consecutive resolutions, a contingency (crosstab) heatmap shows how
parent clusters split into child clusters as resolution increases. Heatmaps and
an Excel workbook are written to `transitions/`.

## CPU vs GPU backend divergence

**This matters for interpretation.** The two backends are *not* numerically
equivalent:

- **CPU** uses `pynndescent` for kNN and `leidenalg`/`igraph` for Leiden
  (`flavor="igraph"`, 2 iterations, undirected).
- **GPU** uses cuML kNN and cuGraph Leiden.

Consequences observed on the pbmc3k benchmark:

- The two backends produce **different cluster counts at identical resolution
  values**.
- GPU Jaccard scores are **systematically lower**, and the gap tends to widen
  with resolution.
- Therefore CPU and GPU results are **not directly comparable**. Treat them as
  distinct backends with distinct result distributions, and keep the backend
  fixed within an analysis.

**Performance.** GPU acceleration is counterproductive on small datasets
(< ~50k cells): the bootstrap consists of many small, sequential
re-clusterings whose host/device transfer and launch overhead dominate. The
bootstrap therefore **auto-routes to CPU** below `--gpu-cell-threshold`
(default 50,000 cells) when the stability backend is `auto`. GPU is expected to
pay off on larger datasets.

## Known limitations

- **GPU reproducibility** — cuGraph Leiden is currently **unseeded**, so GPU
  bootstrap results vary run-to-run. Seeding is on the roadmap; for reproducible
  runs today, use the CPU backend.
- **Backend comparability** — see above; do not mix CPU and GPU cluster labels
  in a single comparison.
- **GPU warm-up** — a pre-loop kernel warm-up reduces, but may not fully remove,
  one-off JIT/allocation cost charged to the first timed resolution.

## Performance metrics

When profiling is enabled (default; disable with `--skip-perf`), per-stage
wall-clock and per-resolution bootstrap throughput are written to
`performance/` as CSV, JSON, and PNG, and embedded in the HTML report header.
