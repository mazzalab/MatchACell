# scParadise annotation

MatchACell runs the **scAdam v2** cell type annotator from
[scParadise](https://github.com/Chechekhins/scParadise) as an independent Step-2
Snakemake rule, alongside the existing annotators. It uses a local pretrained
or custom model and the clustered `.h5ad` and MatchA Verdict from Step 1.
The rule is included in `annotation` only when
`matchacell_annotation.scparadise.model_dir` is nonempty. Existing configurations
without a `scparadise` section continue to skip it.

## Model and environment

The isolated rule environment, `workflow/envs/scparadise.yaml`, pins
`scparadise==1.1.0`. Snakemake creates it automatically with conda enabled.
For model preparation or standalone execution, create it explicitly:

```bash
conda env create -f workflow/envs/scparadise.yaml
conda activate matchacell-scparadise
python - <<'PY'
from scparadise import scadam
print(scadam.available_models())
# Example for human PBMC data; select a model matching your tissue/species.
scadam.download_model('Human_PBMC', save_path='models')
PY
```

Use the resulting `models/Human_PBMC_scAdam` directory, containing
`model_v2.pth`. Custom models must use the same v2 format and a folder name
ending in `scAdam`, which scParadise's prediction API expects. Legacy TabNet v1
models are not supported by this integration. Unknown detection additionally
requires the model's `unknown_detector.json`.

The model is shared across samples and is an explicit workflow input; annotation
does not download or train models. Keep downloaded checkpoints outside Git.

## Configure and run

Set these entries in your working configuration (alongside the other annotators):

```yaml
matchacell_annotation:
  scparadise:
    model_dir: "models/Human_PBMC_scAdam"
    layer: "lognorm"
    label_level: 0
    threshold: 0.0
    detect_unknown: false
    unknown_method: "voting"
    unknown_threshold: 2.0
    batch_size: 256
    device: "auto"
    n_threads: 8
    seed: 0
```

From the usual MatchACell launcher environment:

```bash
python run.py -w annotation -c config.yaml -q 8 -n
python run.py -w annotation -c config.yaml -q 8
```

To run just this annotator with Snakemake, target its output file using the
configured output directory and sample ID:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 8 --use-conda \
  /path/to/output/results/pbmc3k/matchacell/annotation/scParadise/scparadise_annotated.h5ad
```

Or run the script in the scParadise environment after Step 1:

```bash
python workflow/scripts/scparadise_tool.py \
  --input /path/to/clustered_multi_resolution.h5ad \
  --verdict /path/to/MatchA_Verdict.txt \
  --model_dir models/Human_PBMC_scAdam \
  --output /path/to/annotation/scParadise --device cpu
```

## Expression, annotation levels and unknown cells

- `layer` selects **log-normalized expression**, defaulting to Step 1's
  `lognorm` layer. Use `""` to select `.X`. This follows the upstream
  [normalization guidance](https://scparadise.readthedocs.io/en/latest/tutorials/notebooks/scAdam/scAdam_predict.html).
  Do not select raw counts or scaled/centered expression. The script does not
  normalize again or change the input's expression, counts, `.raw`, or embeddings.
- Model and input gene identifiers must match, including species and naming
  convention. No overlap is an error; less than 80% emits a warning. Missing
  model genes are zero-filled by scAdam. Step 1's default HVG selection can
  discard model genes: set `matchacell_cluster_stability.extra: "--use-hvg no"`
  and regenerate Step 1 when you need to retain them. Targeted spatial panels
  can also have low overlap; interpret their predictions accordingly.
- `label_level: 0` selects the finest available level for reporting; positive
  integers select a particular level. All levels and probabilities are retained.
- `threshold` is a probability cutoff for the selected label in `cell_type_pred`.
  Values below it become `Unknown`; `0.0` disables this additional filter.
- `detect_unknown` enables scAdam's own detector. `unknown_method` accepts
  `voting`, `gradient`, `entropy`, `distance`, or `combined`.
  `unknown_threshold` is its **Ashman separation threshold**, not a probability.
  The other detector parameters use scParadise defaults. See the upstream
  [prediction API](https://scparadise.readthedocs.io/en/latest/api/generated/scparadise.scadam.predict.html).
- `device` accepts `auto`, `cpu`, or `cuda`. `auto` uses an available CUDA device;
  `cuda` fails clearly if none is available. GPU use requires a CUDA-capable
  PyTorch installation and driver in the rule environment. Snakemake caps CPU
  threads to its allocated cores. scAdam densifies shared model genes internally;
  `batch_size` limits inference batches, not the complete expression allocation.

## Outputs

Files are written under
`<output_dir>/results/<sample>/matchacell/annotation/scParadise/`:

- `scparadise_annotated.h5ad`: the full clustered object with
  `scparadise_celltype_l1`, `scparadise_celltype_l2`, etc., their `_probability`
  columns, and `predicted_celltype`, `prediction_probability`, `cell_type_pred`
  for the selected level. Unknown flags/scores, when returned, have the
  `scparadise_` prefix. `.uns['scparadise']` records model path, version,
  parameters, selected level, and model gene overlap.
- `cluster_to_celltype_mapping.xlsx`: counts of selected cell labels by Verdict
  cluster, plus `dominant_cell_type` and `dominant_pct`, matching scANVI's layout.
- UMAP plots of the selected labels and a panel per label, using the shared
  plotting helper. If Step 1 has no UMAP, annotation and the table still succeed.

The tests cover input handling, label selection, unknown labels, and preservation
of the clustered object. With scParadise installed, they also run real CPU
inference on a tiny synthetic checkpoint; this validates integration, not
biological accuracy on a reference atlas.
