# CyteType

Sends each Verdict cluster's marker genes to the
[CyteType](https://cytetype.nygen.io) AI annotation service and writes back its
cluster annotations, including Cell Ontology terms.

| | |
| --- | --- |
| Runs | when a CyteType API token is set |
| Labels | per cluster |
| Rule / script | `workflow/rules/cytetype.smk` → `workflow/scripts/cytetype_tool.py` |
| Environment | `workflow/envs/cytetype.yaml` (Python 3.12, `cytetype` from PyPI) |
| Needs | a CyteType account and API token; internet access where the rule runs |

## Configuration

```yaml
matchacell_annotation:
  cytetype:
    api_token: ""      # keep empty; use the CYTETYPE_API_TOKEN environment variable
    study_context: "Human PBMC scRNA-seq data (10x Genomics), ..."
    title: "MatchACell run"
    run_label: "v1"
    n_top_genes: 100   # marker genes per cluster sent to CyteType
```

`config.yaml` is committed to git, so leave `api_token` empty and export the
token instead:

```bash
export CYTETYPE_API_TOKEN=cyt_...
./run.py -w annotation -c config.yaml -q 8
```

The workflow schedules CyteType only when a token is set in `config.yaml` or in
that variable. The script can also use credentials saved by `cytetype setup`,
but those alone don't enable the rule.

`study_context` is free text describing the samples, such as tissue, species and
technology; CyteType uses it to interpret the markers.

## How it works

1. Ranks marker genes for each Verdict cluster with a Wilcoxon test on the
   log-normalized `.X`, stored as `rank_genes_<leiden_res>`.
2. Builds a query from each cluster's top `n_top_genes` markers, the study
   context and the run metadata, saves it as `query.json`, and submits it to
   CyteType.
3. Writes the returned annotation for each cluster onto its cells.

## Outputs

In `annotation/CyteType/`:

- `cytetype_annotated.h5ad`: `cytetype_annotation_<leiden_res>`, plus
  `cytetype_cellOntologyTerm_<leiden_res>`,
  `cytetype_cellOntologyTermID_<leiden_res>` and
  `cytetype_cellState_<leiden_res>` when CyteType returns them.
- `cluster_to_celltype_mapping.xlsx`: those columns, one row per cluster.
- `query.json`: the query submitted to CyteType.
- UMAPs of the annotation (overall and one panel per label).

## Notes

- Your data leaves the machine. Check `query.json` against your data-sharing
  rules before running on sensitive samples.
