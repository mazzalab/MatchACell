# reference_pbmc.h5ad

Annotated single-cell reference used to test ScInfeR's reference-mapping
annotation (`matchacell_annotation.scinfer.reference_file` in
[config.yaml](../config.yaml)).

## What it is

The classic 10x Genomics PBMC 3k dataset, pre-processed and annotated via
the standard Scanpy/Seurat clustering tutorial: 2,638 peripheral blood
mononuclear cells (after QC filtering from the original ~2,700), 1,838
genes, log-normalized expression in `.X`.

| | |
| --- | --- |
| `obs['louvain']` | cell type label (8 classes, see below) |
| `obsm['X_umap']` | 2D UMAP embedding, required by `ScInfeR::get_marker_from_ref_matrix()` |
| `obsm['X_pca']`, `X_tsne`, `X_draw_graph_fr` | other embeddings from the tutorial, unused by ScInfeR |
| `obsp['distances']`, `connectivities'` | neighbor graph from the tutorial |

Cell type counts in `louvain`:

```
CD4 T cells          1144
CD14+ Monocytes       480
B cells               342
CD8 T cells           316
NK cells              154
FCGR3A+ Monocytes     150
Dendritic cells        37
Megakaryocytes         15
```

## Why this dataset

MatchACell's own test query (`data/pbmc3k.h5ad`, `data/pbmc3k_raw.h5ad`) has
no cell-type labels and no UMAP — neither `ScInfeR.R`'s `--reference`
requirement is met, so it can't serve as its own reference. The 10x
Datasets portal (e.g. the 5k Human Donor2 PBMC listing) only ships raw/
filtered count matrices, not cell-type labels, so it doesn't help either.
`pbmc3k_processed` is the smallest dataset that already carries real
cell-type labels and a UMAP for the same PBMC cell cohort used elsewhere in
the test config, with no manual annotation step required.

## How it was obtained

```python
import scanpy as sc
a = sc.datasets.pbmc3k_processed()   # scanpy's built-in copy of the tutorial dataset
a.write_h5ad("data/reference_pbmc.h5ad")
```

`sc.datasets.pbmc3k_processed()` also writes its own cache copy to
`data/pbmc3k_processed.h5ad` (relative to the working directory at call
time); that file was deleted after confirming it was byte-for-byte
equivalent in content to `reference_pbmc.h5ad` (same `.obs`, `.var`, `.X`,
and `obsm['X_umap']`), to avoid keeping two copies of the same data.

## Using it

`config.yaml`:

```yaml
matchacell_annotation:
  scinfer:
    reference_file: "data/reference_pbmc.h5ad"
    ref_label_col: "louvain"   # not "cell_type" -- that's the column name in this dataset
```

The `scinfer` Snakemake rule converts this `.h5ad` to a Seurat `.rds` via
`workflow/scripts/h52rds.R` before handing it to `ScInfeR.R`.
