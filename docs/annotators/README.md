# Step 2 — Annotators

The `annotation` target labels the clustering chosen in Step 1 with up to eight
annotators. Each one is its own Snakemake rule, with its own script and conda
environment. They don't depend on one another, so Snakemake runs them in
parallel up to `--cores`.

```bash
./run.py -w annotation -c config.yaml -q 16
```

## The annotators

| Annotator | Approach | Labels | Runs when | Needs |
| --- | --- | --- | --- | --- |
| [ScoreGenes](score-genes.md) | Marker-gene signatures scored with Scanpy `score_genes` | per cluster | always | the signature workbook |
| [AddModuleScore](addmodulescore.md) | Marker-gene signatures scored with Seurat `AddModuleScore` (R) | per cluster | always | the signature workbook |
| [CIA](cia.md) | Marker-gene signatures, each cell classified by CIA | per cell, summarised per cluster | always | the signature workbook |
| [CellTypist](celltypist.md) | Pretrained CellTypist models | per cell, majority vote per cluster | always | internet access on the first run |
| [scANVI](scanvi.md) | Reference mapping with scvi-tools | per cell | `scanvi.reference_file` is set | an annotated reference `.h5ad` |
| [scParadise](scparadise.md) | Pretrained scAdam model | per cell | `scparadise.model_dir` is set | a downloaded scAdam model |
| [CyteType](cytetype.md) | CyteType AI annotation service | per cluster | a CyteType API token is set | a CyteType account, internet access |
| [CellTypeAI](celltypeai.md) | Local LLM ensemble served by Ollama | per cluster | `celltypeai.tissue` is not empty | a running Ollama server |

Each annotator speaks its own vocabulary: the sheet names of your workbook, a
model's or reference's labels, or free text. Reconciling them is the job of the
planned Step 3.

## What every annotator shares

**Input.** Each rule reads Step 1's `clustered_multi_resolution.h5ad`, which
holds log-normalized expression in `.X` plus `counts` and `lognorm` layers.
AddModuleScore reads its Seurat companion, `clustered_multi_resolution.rds`,
instead.

**Resolution.** Each script reads `MatchA_Verdict.txt` and annotates the
`leiden_<res>` column the Verdict recommends. Rerunning Step 1 with different
settings can change that resolution, and with it every annotation.

**Output folder.** Everything goes to
`<output_dir>/results/<sample>/matchacell/annotation/<Annotator>/`: the
annotated `.h5ad` (the rule's tracked output), a per-cluster Excel table, and
UMAP plots. Where the labels end up:

| Annotator | Label column in the `.h5ad` | Cluster table |
| --- | --- | --- |
| ScoreGenes | `cell_type_pred`, plus one score column per signature | `cluster_annotation_summary.xlsx` |
| AddModuleScore | `cell_type_pred`, plus one score column per signature | `cluster_annotation_summary.xlsx` |
| CIA | `CIA prediction default` | `cluster_to_celltype_mapping.xlsx` |
| CellTypist | `celltypist_label_majority_voting_<model>` | `cluster_to_celltype_mapping.xlsx` |
| scANVI | `cell_type_pred`, `prediction_probability` | `cluster_to_celltype_mapping.xlsx` |
| scParadise | `cell_type_pred`, `prediction_probability` | `cluster_to_celltype_mapping.xlsx` |
| CyteType | `cytetype_annotation_<leiden_res>` | `cluster_to_celltype_mapping.xlsx` |
| CellTypeAI | `cell_type_ai` | `cluster_to_celltype_mapping.xlsx` |

ScoreGenes, AddModuleScore, scANVI and scParadise write `Unknown` when nothing
passes their threshold. CIA writes `Unassigned` for cells it can't assign and
`Unknown` for clusters where no label reaches its threshold.

**Enabling annotators.** ScoreGenes, AddModuleScore, CIA and CellTypist always
run with `annotation`; the other four are skipped while their prerequisite in
the table above is empty. Keep every annotator's section in `config.yaml`, plus
`annot_file`, even for results you ignore: the rules of ScoreGenes, CIA,
CellTypist, AddModuleScore, CyteType and CellTypeAI read their sections when the
workflow is parsed, and four rules list the workbook as an input, so the file
must exist.

## Keep every gene for annotation

With Step 1's default `--use-hvg auto`, a matrix with more than 2,000 genes is
reduced to its 2,000 most variable genes, and that reduced object is what every
annotator reads. Pretrained models and marker signatures then lose most of their
genes: on a PBMC sample, the T-cell marker *CD3D* was among those dropped. Before
running Step 2, keep all genes and rerun Step 1:

```yaml
matchacell_cluster_stability:
  extra: "--use-hvg no"
```

PCA and clustering then use all genes too. This also avoids an
[AddModuleScore failure on small gene sets](addmodulescore.md#notes).

## The signature workbook

ScoreGenes, AddModuleScore and CIA read the same Excel workbook,
`matchacell_annotation.annot_file`:

- one sheet per cell type or state; **the sheet name is the label** those
  annotators report;
- gene symbols in the first column, one per row, **no header row**;
- symbols must match the dataset's gene names, including species and
  capitalisation.

Genes missing from the dataset are dropped, and so are sheets left empty. CIA
also ignores signatures with fewer than 3 genes. Excel limits sheet names to 31
characters and doesn't allow `: \ / ? * [ ]`.

For example, a workbook with three sheets:

| Sheet `CD14 Mono` | Sheet `NK` | Sheet `B cells` |
| --- | --- | --- |
| CD14 | NCAM1 | MS4A1 |
| LYZ | KLRF1 | CD79A |
| S100A8 | GNLY | CD79B |
| FCN1 | NKG7 | BANK1 |

## Running a single annotator

`run.py` targets whole steps. To run one annotator, ask Snakemake for its output
file; it runs Step 1 first if needed:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 8 --use-conda \
  /path/to/output/results/pbmc3k/matchacell/annotation/CellTypist/celltypist_annotated.h5ad
```

To add an annotator, see [Adding an annotator](../pipeline.md#adding-an-annotator).
