# AddModuleScore

Scores the same marker-gene signatures as [ScoreGenes](score-genes.md), but with
Seurat's [`AddModuleScore`](https://satijalab.org/seurat/reference/addmodulescore)
in R. The two implementations of this control-gene scoring are kept separate so
they can confirm each other instead of sharing one implementation's blind spots.

| | |
| --- | --- |
| Runs | always, as part of `annotation` |
| Labels | per cluster |
| Rule / script | `workflow/rules/addmodulescore.smk` → `workflow/scripts/addmodulescore_tool.R`, then `rds2h5.R` |
| Environment | `workflow/envs/addmodulescore.yaml` (R 4.3, Seurat ≥ 5, zellkonverter) |
| Needs | the [signature workbook](README.md#the-signature-workbook) |

## Configuration

```yaml
matchacell_annotation:
  annot_file: "signatures.xlsx"
  addmodulescore:
    threshold: 0.1   # minimum cluster-mean module score for a signature to be called
```

## How it works

1. Reads Step 1's Seurat companion, `clustered_multi_resolution.rds` (built by
   the `matchacell_cluster_stability_rds` rule), and sets the Verdict
   resolution as the cell identities.
2. Keeps each signature's genes that are present in the object.
3. Runs `AddModuleScore` once per signature with 24 expression bins. It asks for
   100 control genes per bin, capped at the average bin size (genes ÷ 24) for
   small gene sets.
4. Averages each module score per cluster, then calls signatures exactly like
   ScoreGenes: every signature reaching `threshold`, strongest first, joined
   with `, `; otherwise `Unknown`.
5. Saves the Seurat object and converts it back to `.h5ad` with `rds2h5.R`.

## Outputs

In `annotation/AddModuleScore/`:

- `addmodulescore_annotated.h5ad`: one module-score column per signature and
  `cell_type_pred`.
- `rds/addmodulescore_annotated.rds`: the same object, for use in R.
- `cluster_annotation_summary.xlsx`: per cluster, the mean score of each
  signature, `final_annotation` and `highest_score`.
- `violin/<signature>_violin.png` and `Heatmap_Signature_Validation.png`.
- `<leiden_res>_leiden_umap.png`, `<signature>_umap.png` for each signature, and
  `<leiden_res>_cell_type_pred_umap.png`.

## Notes

- **Small gene sets can fail.** On a 2,000-gene object (Step 1's default of
  keeping only highly variable genes), Seurat stopped with
  `cannot take a sample larger than the population when 'replace = FALSE'`:
  an expression bin held fewer genes than the control genes requested. Keep all
  genes, as described in
  [Keep every gene for annotation](README.md#keep-every-gene-for-annotation).
- The UMAP plots use Step 1's `X_umap` reduction. If it isn't found, the script
  warns and falls back to the object's first reduction.
