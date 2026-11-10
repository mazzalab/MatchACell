# Step 3: consensus across annotators

Step 3 reconciles the annotators of Step 2 on the Cell Ontology (CL). It works
on the clusters of the MatchA Verdict resolution and writes one consensus label
per cluster, with how much of the vote backs it. It ends with
`MatchACell_full_report.html`, one report for Steps 1 to 3.

```bash
./run.py -w consensus -c config.yaml -q 16   # runs Steps 1 and 2 first if needed
```

## How the vote works

The vote is the ontology vote of popV (Ergen et al., *Nat Genet* 2024), applied
to the clusters of Step 1.

1. **Each annotator's call per cell.** Cells labelled `Unknown`, `Unassigned` or
   left empty abstain.
   - ScoreGenes and AddModuleScore label a cluster with every signature above
     their threshold, best score first (`CD14 Mono, FCGR3A Mono`): the call is
     the first one.
   - CellTypist runs several models. Per cell, when the models' terms lie on one
     lineage (Low: "classical monocyte", High: "monocyte") the call is the most
     specific one; when they conflict (Low: "NK cell", High: "T cell") the cell's
     vote is split between them.
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
3. **Votes per cluster.** In a cluster, an annotator gives each term the share
   of the cluster's cells that carry it. Abstaining cells still count in the
   cluster size, so an annotator that says `Unknown` for most of a cluster
   votes with what is left.
4. **Weights.** With `vote_by: family` (the default) each family of annotators
   has one vote, split equally between its members that labelled the cluster:

   | Family | Annotators |
   | --- | --- |
   | marker signatures | ScoreGenes, AddModuleScore, CIA |
   | reference models | CellTypist, scANVI, scParadise |
   | LLM | CyteType, CellTypeAI |

   Annotators of one family share their evidence: the three signature methods
   score the same workbook, and ScoreGenes and AddModuleScore compute nearly the
   same score. On pbmc3k they gave the same label to 100% of the cells. With one
   vote each they would always outvote a single reference model.
   `vote_by: annotator` gives every annotator one vote. `weights` multiplies an
   annotator's vote in either mode.
5. **Hierarchy.** A vote for a term also supports every `is_a` ancestor of it.
   "CD4-positive, alpha-beta T cell" and "CD8-positive, alpha-beta T cell"
   disagree, but both support "T cell".
6. **Consensus.** The most specific term supported by more than `min_agreement`
   of the weighted vote (default: more than half).

| `status` | Meaning |
| --- | --- |
| `agreed` | The consensus term is a label some annotator gave. |
| `partial` | Annotators only agree on a coarser common ancestor (e.g. CD4 vs CD8 → T cell). |
| `none` | No term reaches the bar, e.g. unrelated free-text labels. |

## Validation against known labels

When the input has a column of known labels, for example from a published
annotation, set `truth_col` to its name. Step 3 then scores, per cell, every
annotator and the consensus. Known labels are mapped to CL with the same label
maps, and each call is:

- **correct**: the known term, or a more specific term under it;
- **coarser**: an ancestor of the known term, right but less specific;
- **wrong**: a term in another lineage;
- **no call**: the annotator abstained, or the cluster has no consensus.

`Cluster ceiling` labels each cluster with its most common known label. No
annotation that gives a whole cluster one label can do better at this
resolution.

On pbmc3k, scored against the labels of the scanpy/Seurat tutorial (2,327
cells after QC and scDblFinder, Verdict leiden_1.5, 11 clusters):

| | correct | coarser | wrong |
| --- | --- | --- | --- |
| ScoreGenes / AddModuleScore | 89% | 0% | 11% |
| CIA | 87% | 0% | 13% |
| CellTypist | 89% | 2% | 9% |
| **Consensus** (either vote mode) | 89% | 0% | 11% |
| Cluster ceiling | 93% | 0% | 7% |

All 11 clusters are `agreed`. Nearly every error is one cluster of 124 CD8 T
cells that ScoreGenes, AddModuleScore and CellTypist's fine model call NK
(these cells express *NKG7*): the consensus can't be right when the annotators
are wrong the same way. The 3 cells left of the tutorial's `Megakaryocytes`
are too few to form a cluster, so only CIA, which labels single cells, gets
them.

At the finer leiden_2.0 (before the Verdict counted cells in unstable
clusters), the two vote modes differed: per family 83% correct, 8% coarser
and 8% wrong; per annotator 89%, 0% and 11%. Pooling by family traded exact
calls for fewer wrong ones by falling back to the term the families share.
The tutorial labels themselves come from marker genes, which favours the
marker methods.

## Label maps

A label map is a TSV with `label` and `cl_id` columns (other columns are
ignored). `matchacell_consensus.label_maps` lists the maps to use; a later map
overrides an earlier one, so put your own last. Every id must exist in the
ontology, or the step stops and names the bad ids.

| Map | Covers |
| --- | --- |
| `workflow/resources/label_maps/celltypist_immune.tsv` | The labels of CellTypist's `Immune_All_Low` and `Immune_All_High` models, with the CL ids of the [CellTypist encyclopedia](https://www.celltypist.org/encyclopedia). Coarse labels without their own id get the most specific common ancestor of their fine labels. Regenerate with `tools/build_celltypist_label_map.py`. |
| `workflow/resources/label_maps/pbmc.tsv` | Common PBMC names (`CD4 T`, `CD14 Mono`, `FCGR3A Mono`, …, and the Seurat/scanpy pbmc3k tutorial names). CD14 and FCGR3A monocytes map to *classical* and *non-classical monocyte*, the terms CellTypist uses: in CL, "CD14-low, CD16-positive monocyte" sits under "CD14-positive monocyte", which would make an FCGR3A vote support the CD14 label. It also maps `Megakaryocytes` and CellTypist's `Megakaryocytes/platelets` to *platelet*: in blood these are platelets, and CL doesn't put platelets under megakaryocytes. |

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
  vote_by: "family"     # or "annotator"
  weights: {}           # e.g. {CellTypist: 2}
  exclude: []           # annotators left out of the vote, e.g. ["CellTypeAI"]
  truth_col: ""         # column of the input .obs with known labels
```

The ontology is the CL release pinned in `workflow/rules/consensus.smk`
(`cl-basic.obo`, the edition with `is_a` relations only), downloaded once into
`<output_dir>/resources/` by a rule that runs on the submitting node. Without
internet there, download it elsewhere and set `ontology`.

## Outputs

In `results/<sample>/matchacell/`:

| File | Content |
| --- | --- |
| `MatchACell_full_report.html` | One self-contained page for Steps 1 to 3 (below). |
| `consensus/consensus_by_cluster.csv` | One row per cluster: `n_cells`, Step 1 `stability`, `consensus_label` (in an annotator's own words when one used that term, else the CL name), `consensus_cl_id`, `consensus_cl_name`, `agreement` (share of the vote that supports the term), `compatibility` (adds the votes for coarser terms, which don't contradict it), `status`, the most voted single term (`top_label`), each annotator's share of cells that support the consensus (`support:<annotator>`) and each label set's leading label (`call:<voter>`). |
| `consensus/annotation_long.csv` | Every label of every annotator per cluster: its CL term, how it was mapped, its share of the cluster and how it relates to the consensus (`same`, `finer`, `coarser`, `conflict`). |
| `consensus/label_vocabulary.csv` | Every label of every annotator, its CL term, how it was mapped and how many cells carry it. |
| `consensus/annotator_agreement.csv` | For every pair of label sets, over the cells both label: the share with the same CL term, and with compatible terms (one an ancestor of the other). |
| `consensus/cell_calls.csv.gz` | Per cell: the cluster, every annotator's label and the consensus. |
| `consensus/consensus_settings.json` | Resolution, ontology release, label maps, vote mode and weights of the run. |
| `consensus/validation.csv`, `validation_by_label.csv` | With `truth_col`: the scores above, overall and per known label. |
| `consensus/consensus_annotated.h5ad` | The Step 1 object with `consensus_label`, `consensus_cl_id`, `consensus_cl_name`, `consensus_agreement` and `consensus_status` in `.obs`. |

## The report

`MatchACell_full_report.html` opens offline (plotly.js is inlined) and follows
the steps:

- **Overview**: cells kept, the Verdict, the annotators, the consensus outcome,
  and an embedding you can colour by any Leiden resolution, annotator, the
  consensus, the agreement, the stability or a QC metric. Labels that map to the
  same CL term get the same colour in every annotator.
- **Step 1**: cells removed per QC criterion and the QC figures; the stability
  of every cluster at every resolution and the Verdict; how clusters split
  between resolutions.
- **Step 2**: what each annotator labelled and how much of it maps to CL; every
  annotator's leading call per cluster, marked as agreeing with the consensus,
  coarser, or conflicting; how often each pair of annotators agrees; the
  signature scores behind ScoreGenes and AddModuleScore.
- **Step 3**: the consensus per cluster; stability against agreement, to see
  which clusters to trust; the marker genes of the signature workbook per
  cluster; the validation when `truth_col` is set.
- **Methods and settings**, with references, and the output files.

## Reading the result

- Low `agreement` on an unstable cluster (low Step 1 `stability`) often means the
  cluster mixes populations; a coarser resolution may annotate better.
- `partial` with a high `agreement` means the annotators agree on the lineage
  but not the subtype: check the calls and the marker genes in the report.
- The consensus can't be better than its annotators: when every family is wrong
  the same way (the CD8 cluster above), it is wrong too. The marker plot is the
  check.

## Not done yet

- A per-cell consensus for annotators that label single cells (popV votes per
  cell; MatchACell votes per cluster).
- Mapping free-text labels (LLM annotators) beyond exact names and synonyms.
