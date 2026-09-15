# ScoreGenes

Scores each marker-gene signature with Scanpy's
[`score_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html)
and labels each Verdict cluster with the signatures it scores highly on.

| | |
| --- | --- |
| Runs | always, as part of `annotation` |
| Labels | per cluster |
| Rule / script | `workflow/rules/score_genes.smk` → `workflow/scripts/score_genes.py` |
| Environment | `workflow/envs/matchacell.yaml` (the Step 1 environment) |
| Needs | the [signature workbook](README.md#the-signature-workbook) |

## Configuration

```yaml
matchacell_annotation:
  annot_file: "signatures.xlsx"
  score_genes:
    threshold: 0.1   # minimum cluster-mean score for a signature to be called
```

## How it works

1. Reads the workbook and keeps each signature's genes that are present in the
   data.
2. Runs `sc.tl.score_genes` once per signature on the log-normalized `.X`. A
   cell's score is the mean expression of the signature genes minus the mean of
   control genes drawn from the same expression bins. Each score is stored in a
   column named after its sheet.
3. Averages every score per Verdict cluster.
4. Calls every signature whose cluster mean reaches `threshold`, strongest
   first, joined with `, ` (for example `T cells, CD4 T`). A cluster where none
   does is `Unknown`.
5. Writes each cluster's call onto its cells as `cell_type_pred`.

## Outputs

In `annotation/ScoreGenes/`:

- `score_genes_annotated.h5ad`: the clustered object with one score column per
  signature and `cell_type_pred`.
- `cluster_annotation_summary.xlsx`: per cluster, the mean score of each
  signature, `final_annotation`, and `highest_score` (the best passing score, or
  the best score overall for `Unknown` clusters).
- `violin/<signature>_violin.png`: the score per cluster, with the global median
  as a dashed line.
- `Heatmap_Signature_Validation.png`: mean score per cluster and signature.
- A dot plot of signature scores per cluster, UMAPs of the clusters and of each
  signature, and UMAPs of `cell_type_pred` (overall and one panel per label).

## Notes

- `threshold` is on the score's own scale, which depends on the signature and
  the dataset. Check the heatmap before relying on the default of 0.1.
- Related types often both pass, such as a general T-cell signature and a CD4 T
  signature, so a cluster can carry several labels.
- [AddModuleScore](addmodulescore.md) computes the same kind of score in R with
  Seurat, as an independent check.
