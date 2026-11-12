# CellWhisperer annotation

MatchACell runs [CellWhisperer](https://github.com/epigen/cellwhisperer)
(Schaefer et al., *Nature Biotechnology* 2025) as an independent Step-2 rule.
CellWhisperer is a CLIP-style model: it embeds transcriptomes (through
Geneformer) and free text into one space. Each cell is scored against every
candidate label, and it gets the label with the highest score.

Both embeddings run on the local model. No data leaves the machine, and the
hosted CellWhisperer API is not used.

The rule is in `annotation` only when `matchacell_annotation.cellwhisperer.repo_dir`
is set. Configurations without a `cellwhisperer` section skip it.

## How it labels cells

- **Candidates.** The labels in `labels`. If `labels` is empty, the sheet names
  of `annot_file` are used, so the vocabulary matches the signature annotators.
- **Text.** Each label is put into `template`, e.g. `"{label} from human pancreas"`.
  Tissue context helps. On a pancreas test set (600 cells, 12 types) it raised
  per-cell accuracy from 0.43 to 0.50.
- **Per cell.** The softmax over the candidates of the CLIP logits.
  `cell_type_pred` is the best label and `prediction_probability` its probability.
  Below `threshold`, the cell is `Unknown`.
- **Per cluster.** `cluster_to_celltype_mapping.xlsx` gives, for each cluster
  of the Verdict resolution:
  - the majority of the cell calls, and its fraction;
  - the call of the cluster's mean embedding (`centroid_label`).

CellWhisperer is zero-shot. It only chooses among the labels you give, and it
can't detect a type that is not in the list. On the pancreas test set, it
separated endocrine, exocrine, endothelial and immune cells. Close subtypes were
weak (alpha vs beta vs delta cells, TFF1+ vs TUBA1A+ ductal cells). Use it as one
vote in the consensus, not on its own.

## Setup (once)

CellWhisperer ships a pixi environment, not a conda one, so this rule runs its
script with `pixi run` in a CellWhisperer checkout. It has no `conda:` directive.

```bash
# 1. pixi (https://pixi.sh), e.g.
brew install pixi            # or: curl -fsSL https://pixi.sh/install.sh | bash

# 2. CellWhisperer with its submodules (UCE is needed even with Geneformer)
git clone --recurse-submodules https://github.com/epigen/cellwhisperer.git
cd cellwhisperer
bash envs/setup_pixi.sh
pixi add openpyxl            # MatchACell reads and writes Excel files

# 3. Model weights (~1.6 GB): the CellWhisperer checkpoint, Geneformer and the
#    gene-symbol map. The simplest way to get them is CellWhisperer's own
#    preprocessing pipeline on any small dataset (see its README), which
#    downloads them into results/models/ and resources/.
```

A checkout cloned without `--recurse-submodules` fails with
`ModuleNotFoundError: No module named 'UCE'`. Then run
`git submodule update --init --recursive` and `pixi reinstall cellwhisperer`.

## Configure and run

```yaml
matchacell_annotation:
  cellwhisperer:
    repo_dir: "/path/to/cellwhisperer"
    model: ""                       # default: <repo_dir>/results/models/jointemb/cellwhisperer_clip_v1.ckpt
    pixi: "pixi"                    # full path if pixi is not on PATH on compute nodes
    labels: []                      # empty -> sheet names of annot_file
    template: "{label} from human pancreas"
    threshold: 0.0
    batch_size: 32
    device: "auto"                  # auto | cpu | cuda
    n_threads: 8
```

```bash
./run.py -w annotation -c config.yaml -q 8
```

To run only this annotator, target its output file:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 8 \
  /path/to/output/results/<sample>/matchacell/annotation/CellWhisperer/cellwhisperer_annotated.h5ad
```

Or run the script directly after Step 1:

```bash
pixi run --manifest-path /path/to/cellwhisperer/pixi.toml \
  python workflow/scripts/cellwhisperer_tool.py \
  --input  /path/to/clustered_multi_resolution.h5ad \
  --verdict /path/to/MatchA_Verdict.txt \
  --output /path/to/annotation/CellWhisperer \
  --model  /path/to/cellwhisperer/results/models/jointemb/cellwhisperer_clip_v1.ckpt \
  --repo_dir /path/to/cellwhisperer \
  --labels "alpha cell" "beta cell" "acinar cell" \
  --template "{label} from human pancreas"
```

## Outputs

In `<output_dir>/results/<sample>/matchacell/annotation/CellWhisperer/`:

| File | Content |
| --- | --- |
| `cellwhisperer_annotated.h5ad` | Step 1 object plus `cell_type_pred`, `prediction_probability`, the cell embeddings in `.obsm["X_cellwhisperer"]`, and the run settings in `.uns["cellwhisperer"]` |
| `cluster_to_celltype_mapping.xlsx` | per-cluster majority and centroid labels |
| `cellwhisperer_cell_scores.csv.gz` | per-cell logit for every candidate label |
| `umap_*.png` | UMAP of the labels, overall and one panel per label |

## Notes

- **Input.** The script embeds the raw `counts` layer. Geneformer normalizes
  and ranks genes in each cell, so keep all genes in Step 1
  (`--use-hvg no`, see [the overview](README.md#keep-every-gene-for-annotation)).
  With fewer than 5,000 genes the script prints a warning.
- **Spatial panels.** Targeted panels (e.g. Xenium, a few hundred genes) give
  Geneformer few genes to rank. Expect weaker labels than on whole-transcriptome
  data.
- **Species.** Human only: the model and the gene mapping are human.
- **Runtime.** On CPU it takes about 1 s per cell (600 cells: ~9 min on an 8-core
  laptop, model loading included). With a GPU (4 GB VRAM is enough) it is much
  faster. The default walltime for `--cluster` is 12 h. Raise it with
  `cluster_resources: {cellwhisperer: {runtime: ...}}` for large samples.
- **Consensus (Step 3).** The label column is `cell_type_pred`.
  `consensus_tool.py` needs a `LABEL_COLUMNS` entry and a `FAMILIES` entry for
  `CellWhisperer`, and `consensus.smk` needs one in `_ANNOTATOR_OUTPUTS`.
