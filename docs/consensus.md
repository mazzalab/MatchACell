# Step 3: consensus across annotators

Step 3 reconciles the annotators of Step 2 on the Cell Ontology (CL). It works
on the clusters of the MatchA Verdict resolution and writes one consensus label
per cluster, with how much of the vote backs it.

```bash
./run.py -w consensus -c config.yaml -q 16   # runs Steps 1 and 2 first if needed
```

## How the vote works

1. **Votes.** Each enabled annotator votes in every cluster with its per-cell
   labels: a label gets the fraction of the cluster's cells that carry it. Cells
   labelled `Unknown`, `Unassigned` or left empty abstain. ScoreGenes and
   AddModuleScore ties such as `CD8 T, CD4 T` split the cell between the labels.
   CellTypist's models share CellTypist's one vote.
2. **Labels to CL terms**, in this order:
   - the CL id the annotator gives itself (CyteType);
   - the [label maps](#label-maps);
   - a label that is itself a CL id (`CL:0000236`);
   - a CL term name or exact synonym, ignoring case, spacing and plurals.
     Abbreviations of up to four characters only match term names: as synonyms
     they are ambiguous (CL lists `DC2` as an old name of plasmacytoid dendritic
     cells).

   Labels that match nothing stay as free text: they only agree with the very
   same label.
3. **Hierarchy.** A vote for a term also supports every `is_a` ancestor of it.
   "CD4-positive, alpha-beta T cell" and "CD8-positive, alpha-beta T cell"
   disagree, but both support "T cell".
4. **Consensus.** The most specific term supported by more than `min_agreement`
   of the annotators that voted in the cluster (default: more than half).

| `status` | Meaning |
| --- | --- |
| `agreed` | The consensus term is a label some annotator gave. |
| `partial` | Annotators only agree on a coarser common ancestor (e.g. CD4 vs CD8 → T cell). |
| `none` | No term reaches the bar, e.g. unrelated free-text labels. |

## Label maps

A label map is a TSV with `label` and `cl_id` columns (other columns are
ignored). `matchacell_consensus.label_maps` lists the maps to use; a later map
overrides an earlier one, so put your own last. Every id must exist in the
ontology, or the step stops and names the bad ids.

| Map | Covers |
| --- | --- |
| `workflow/resources/label_maps/celltypist_immune.tsv` | The labels of CellTypist's `Immune_All_Low` and `Immune_All_High` models, with the CL ids of the [CellTypist encyclopedia](https://www.celltypist.org/encyclopedia). Coarse labels without their own id get the most specific common ancestor of their fine labels. Regenerate with `tools/build_celltypist_label_map.py`. |
| `workflow/resources/label_maps/pbmc.tsv` | Common PBMC names (`CD4 T`, `CD14 Mono`, `FCGR3A Mono`, …, and the Seurat pbmc3k tutorial names). CD14 and FCGR3A monocytes map to *classical* and *non-classical monocyte*, the terms CellTypist uses: in CL, "CD14-low, CD16-positive monocyte" sits under "CD14-positive monocyte", which would make an FCGR3A vote support the CD14 label. It also maps CellTypist's `Megakaryocytes/platelets` to *platelet*: in blood these are platelets, and CL doesn't put platelets under megakaryocytes. |

For your own signature workbook, write one row per sheet name. After a run,
`consensus/label_vocabulary.csv` lists every label each annotator used and
what it was mapped to; rows with `mapped_by = unmapped` are the ones to add.

## Configuration

```yaml
matchacell_consensus:
  ontology: ""          # local cl-basic.obo; empty -> download the pinned release
  label_maps:
    - "workflow/resources/label_maps/celltypist_immune.tsv"
    - "workflow/resources/label_maps/pbmc.tsv"
  min_agreement: 0.5
  exclude: []           # annotators left out of the vote, e.g. ["CellTypeAI"]
```

The ontology is the CL release pinned in `workflow/rules/consensus.smk`
(`cl-basic.obo`, the edition with `is_a` relations only), downloaded once into
`<output_dir>/resources/` by a rule that runs on the submitting node. Without
internet there, download it elsewhere and set `ontology`.

## Outputs

In `results/<sample>/matchacell/consensus/`:

| File | Content |
| --- | --- |
| `consensus_by_cluster.csv` | One row per cluster: `n_cells`, Step 1 `stability`, `consensus_label` (in an annotator's own words when one used that term, else the CL name), `consensus_cl_id`, `consensus_cl_name`, `agreement` (share of the voting annotators that support the term), `compatibility` (share that support it or a coarser term), `status`, the most voted single term (`top_label`) and each annotator's leading label as `call:<annotator>`. |
| `annotation_long.csv` | Every vote: annotator, voter, cluster, label, CL term, how it was mapped, share of the cluster's cells. |
| `label_vocabulary.csv` | Every label of every annotator, its CL term, how it was mapped and how many cells carry it. |
| `consensus_annotated.h5ad` | The Step 1 object with `consensus_label`, `consensus_cl_id`, `consensus_cl_name`, `consensus_agreement` and `consensus_status` in `.obs`. |

## Reading the result

- Low `agreement` on an unstable cluster (low Step 1 `stability`) often means the
  cluster mixes populations; a coarser resolution may annotate better.
- `partial` with a high `agreement` means the annotators agree on the lineage
  but not the subtype: check the `call:` columns and the markers.
- Annotators are weighted equally. An annotator that says `Unknown` for most of
  a cluster votes with what is left.

## Not done yet

- Per-annotator weights (e.g. trusting reference-based methods more).
- Plots and an HTML report.
- Mapping free-text labels (LLM annotators) beyond exact names and synonyms.
- Validation against known labels.
