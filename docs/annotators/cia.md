# CIA

Classifies each cell against the marker-gene signatures with
[CIA](https://github.com/ingmbioinfo/cia) (Cluster Independent Annotation), then
summarises the cell labels per Verdict cluster.

| | |
| --- | --- |
| Runs | always, as part of `annotation` |
| Labels | per cell, summarised per cluster |
| Rule / script | `workflow/rules/cia.smk` → `workflow/scripts/cia_tool.py` |
| Environment | `workflow/envs/cia.yaml` (`cia-python` from PyPI) |
| Needs | the [signature workbook](README.md#the-signature-workbook) |

## Configuration

```yaml
matchacell_annotation:
  annot_file: "signatures.xlsx"
  cia:
    ncpus: 5        # CPUs for CIA's classification
    threshold: 70   # percent of a cluster's cells a label needs to name the cluster
```

## How it works

1. Loads the signatures, prints their pairwise Jaccard similarity (heavily
   overlapping signatures are hard to tell apart), and drops signatures with
   fewer than 3 genes.
2. Runs CIA's `CIA_classify`, which labels each cell with one signature or
   `Unassigned`, in the column `CIA prediction default`.
3. Counts the labels in each Verdict cluster. Every label other than
   `Unassigned` that covers at least `threshold` percent of the cluster's cells
   becomes the cluster's annotation, joined with `, `. A cluster where none does
   is `Unknown`.

## Outputs

In `annotation/CIA/`:

- `cia_annotated.h5ad`: the clustered object with `CIA prediction default` per
  cell.
- `signature_scores_cell.csv`: the per-cell CIA label (despite the file name, it
  holds labels, not scores).
- `cluster_to_celltype_mapping.xlsx`: cell counts per label per cluster, plus
  `CIA_cluster_annotation`.
- UMAPs of the cell labels: overall, with `Unassigned` in grey, and one panel
  per label.

## Notes

- Cells are labelled independently of the clustering, so the cluster table also
  shows how pure each cluster is. With `threshold: 70`, mixed clusters stay
  `Unknown`.
- **Known issue:** the script saves `CIA_annotated.h5ad`, but the rule expects
  `cia_annotated.h5ad`. On case-sensitive filesystems such as Linux, the rule
  then fails with a missing output even though CIA finished. macOS's default
  filesystem hides the mismatch.
