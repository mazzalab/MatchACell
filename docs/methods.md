# Methods & scientific notes

This document describes what the Step-1 engine
(`workflow/scripts/matchacell_cluster_stability.py`) computes and how to read
its output, including the important caveats around the two compute backends.
The Step 2 annotators each describe their method on their own page in
[`annotators/`](annotators/README.md).

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
- **Doublets** (optional, `doublets.method: "scDblFinder"` in the config). The
  `doublets` rule runs scDblFinder (Germain et al. 2021) on the raw counts of
  every cell before QC. It simulates artificial doublets from pairs of cells and
  trains a gradient-boosted classifier to tell them from the real cells; cells
  it calls `doublet` are removed with the other QC criteria (or only flagged,
  `action: "flag"`). In independent benchmarks it is among the most accurate
  doublet detectors (Xi & Li 2021; Germain et al. 2021). Doublets form within
  one capture, so for a file pooling several 10x lanes set
  `doublets.samples_col` to the column naming them. scDblFinder is meant for
  droplet single-cell data, not for imaging-based spatial data such as Xenium.

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
(`--use-hvg auto`). The highly variable genes only drive PCA: the saved object
keeps every gene for the annotators. Leiden is then run across the resolution sweep
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

The result is a per-cluster stability score in `[0, 1]`. This is Hennig's
(2007) cluster-wise bootstrap. His guide for the mean Jaccard (in the
documentation of `fpc::clusterboot`): clusters below 0.6 should not be trusted,
0.6 to 0.75 indicate patterns but not which cells belong to them, 0.75 or more
are valid and stable, 0.85 or more highly stable. MatchACell's 0.60 risk and
0.85 high thresholds and the 0.75 median bar of the Verdict follow it. The PCA is computed once on all cells and kept: the score measures how
much a cluster depends on which cells were sampled, not on gene selection or
the PCA. A boxplot across
resolutions and a per-resolution Excel summary are written to `stability/`.

On CPU every resolution is scored on the same subsamples. Each iteration builds
the subsample's graph once and runs Leiden at every resolution on it; the
subsample's kNN is read off one full-data kNN with `3 × --n-neighbors`
neighbours (a cell's nearest neighbours inside the subsample are its first
full-data neighbours that were drawn). For a given subsample the graph and the
partitions are the same as rebuilding them from scratch with scanpy. Iterations
run on `--n-jobs` worker processes (the rule passes its `threads`), and the
running means are written to `stability/csv/` after each batch, so a job that is
killed keeps what it computed. The workers rebuild each graph with numpy and
igraph only, so they start in seconds. On pbmc3k (2,419 cells, 100 iterations,
8 resolutions) the bootstrap took 14 s on 8 workers, against 147 s for one
resolution at a time, with the same scores.

On GPU the full PCA matrix is staged to device memory once and every iteration
gathers its subsample on-device, avoiding thousands of host→device copies.

## 4. The MatchA Verdict

The Verdict recommends a resolution to carry into annotation.

A resolution is **acceptable** when

- its **median** stability is ≥ `--verdict-min-median` (default 0.75),
- at most `--verdict-max-risk-frac` (default 25%) of its clusters are at risk
  (< `--verdict-risk`, default 0.60), and
- at most `--verdict-max-risk-cells` (default 5%) of its **cells** are in
  clusters at risk.

Among acceptable resolutions, the Verdict picks the **finest** (most clusters):
finer granularity gives more cell-type/state resolution, provided nearly every
cell sits in a reproducible cluster. If none qualify, it falls back to the
resolution with the **highest median stability**.

The share of cells is what decides how much of the data gets labels from
clusters that don't reproduce. Counting clusters alone treats a 12-cell rare
population that is unstable under subsampling like an unstable split of a large
one; it let the Verdict pick resolution 2.0 on pbmc3k (12.7% of the cells in
three unstable CD4 T clusters) and on PBMC 10k (18.1%). With the cell bar,
pbmc3k gets 1.5 (11 clusters, 3.9% of the cells at risk) and PBMC 10k gets 1.0
(20 clusters, median 0.95, 0.2% at risk).

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

- **Doublet detection is off by default.** Doublets often form small clusters
  between two populations; turn on `doublets.method: "scDblFinder"` for droplet
  data. Ambient RNA is not corrected.
- **Stable is not the same as right.** A cluster that merges two cell types
  with similar PCA coordinates can be perfectly stable. Step 3's marker plot and
  the calls of the annotators are the check.

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

## References

- Germain P-L, Lun A, Garcia Meixide C, Macnair W, Robinson MD (2021). Doublet
  identification in single-cell sequencing data using scDblFinder.
  *F1000Research* 10:979.
- Xi NM, Li JJ (2021). Benchmarking computational doublet-detection methods for
  single-cell RNA sequencing data. *Cell Systems* 12:176–194.

- Germain P-L, Sonrel A, Robinson MD (2020). pipeComp, a general framework for
  the evaluation of computational pipelines, reveals performant single cell
  RNA-seq preprocessing tools. *Genome Biology* 21:227. (MAD-based QC)
- Heumos L, Schaar AC, et al. (2023). Best practices for single-cell analysis
  across modalities. *Nature Reviews Genetics* 24:550–572.
- Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing
  well-connected communities. *Scientific Reports* 9:5233.
- Hennig C (2007). Cluster-wise assessment of cluster stability. *Computational
  Statistics & Data Analysis* 52:258–271.
- Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene
  expression data analysis. *Genome Biology* 19:15.
- Ergen C, et al. (2024). Consensus prediction of cell type labels in
  single-cell data with popV. *Nature Genetics* 56:2731–2738. (Step 3's
  ontology vote)
- Diehl AD, et al. (2016). The Cell Ontology 2016: enhanced content,
  modularization, and ontology interoperability. *Journal of Biomedical
  Semantics* 7:44.
