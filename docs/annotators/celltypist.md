# CellTypist

Labels cells with pretrained [CellTypist](https://www.celltypist.org) models and
refines the calls by majority voting within each Verdict cluster.

| | |
| --- | --- |
| Runs | always, as part of `annotation` |
| Labels | per cell, majority vote per cluster |
| Rule / script | `workflow/rules/celltypist.smk` → `workflow/scripts/celltypist_tool.py` |
| Environment | `workflow/envs/celltypist.yaml` (`celltypist` from PyPI) |
| Needs | internet access on the first run; models matching your tissue and species |

## Configuration

```yaml
matchacell_annotation:
  celltypist:
    models: ['Immune_All_Low.pkl', 'Immune_All_High.pkl']
```

List any models from the [CellTypist model collection](https://www.celltypist.org/models);
each runs independently. The defaults are human immune models:
`Immune_All_High` (32 broad types) and `Immune_All_Low` (98 finer types).

## How it works

1. Downloads the CellTypist model collection into CellTypist's model folder
   (`~/.celltypist` by default), skipping files already there. The full
   collection was 62 files, about 120 MB, when last checked.
2. For each model, runs `celltypist.annotate` with majority voting over the
   Verdict clusters: cells are predicted individually, then each cluster takes
   its most common label.
3. Stores each model's majority-voting label as
   `celltypist_label_majority_voting_<model>`.

## Outputs

In `annotation/CellTypist/`:

- `celltypist_annotated.h5ad`: one majority-voting label column per model.
- `cluster_to_celltype_mapping.xlsx`: one row per cluster, one column per model.
- UMAPs of each model's labels (overall and one panel per label).

## Notes

- Only the majority-voting labels are kept; the per-cell predictions and
  CellTypist's confidence scores aren't saved.
- Models need their genes. `Immune_All_Low` and `Immune_All_High` use 6,639
  genes, most of which are gone if Step 1 keeps only highly variable genes; see
  [Keep every gene for annotation](README.md#keep-every-gene-for-annotation).
- The rule lists the signature workbook (`annot_file`) as an input even though
  CellTypist doesn't use it, so that file must exist.
- If compute nodes have no internet access, download the models once
  beforehand from the CellTypist environment:
  `python -c "from celltypist import models; models.download_models()"`.
