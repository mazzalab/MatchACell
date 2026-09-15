# scANVI

Transfers labels from an annotated reference atlas to your cells with
[scANVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html)
from scvi-tools: it trains on the reference, maps your cells onto it, and
predicts a label and a probability for each cell.

| | |
| --- | --- |
| Runs | when `scanvi.reference_file` is set |
| Labels | per cell |
| Rule / script | `workflow/rules/scanvi.smk` → `workflow/scripts/scanvi_tool.py` |
| Environment | `workflow/envs/scanvi.yaml` (scvi-tools, PyTorch, JAX) |
| Needs | an annotated reference `.h5ad` with the same species and gene naming |

## Configuration

```yaml
matchacell_annotation:
  scanvi:
    reference_file: ""               # empty: skip scANVI
    label_col: "free_annotation"     # reference .obs column holding the labels
    unlabeled_category: "Unknown"    # label SCANVI treats as "no label"
    threshold: 0.7                   # below this probability, cell_type_pred is Unknown
    max_epochs_ref: 20               # scVI training on the reference
    max_epochs_scanvi: 20            # SCANVI classifier training
    max_epochs_query: 10             # mapping your cells onto the model
    batch_size: 2048
    n_threads: 10
```

## The reference

- Raw counts, in `.X`, or in `.raw`, which is used when present.
- A label column, `label_col`. Missing labels, including the literal string
  `nan` that some atlases store, become `unlabeled_category`.
- The same species and gene symbols as your data. Only shared genes are used,
  and no shared genes is an error.

## How it works

1. Uses raw counts for both your cells (the `counts` layer) and the reference.
2. Keeps the genes the two have in common.
3. Trains scVI on the reference: 2 layers, 30 latent dimensions,
   negative-binomial likelihood, early stopping, no batch key.
4. Converts it into a SCANVI classifier and trains it on the reference labels.
5. Maps your cells onto the model with `load_query_data` and trains briefly.
6. Predicts `predicted_celltype` and `prediction_probability` (the highest class
   probability) for each cell, and writes `cell_type_pred`: the prediction, or
   `Unknown` below `threshold`.

## Outputs

In `annotation/scANVI/`:

- `scanvi_annotated.h5ad`: `predicted_celltype`, `prediction_probability` and
  `cell_type_pred` per cell; `.X` stays log-normalized.
- `cluster_to_celltype_mapping.xlsx`: cell counts per label per cluster, with
  `dominant_cell_type` and `dominant_pct`.
- UMAPs of `cell_type_pred` (overall and one panel per label).

## Notes

- Batches within the reference (donors, studies) aren't modelled, because no
  batch key is set.
- Training time grows with the size of the reference and the number of epochs.
- The script's own default `threshold` is 0.5; through the workflow, the config
  value (0.7 by default) applies.
