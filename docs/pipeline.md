# The pipeline

This document covers the Snakemake side: the launcher, the configuration, the
rules behind each step, and how to add an annotator.

## Architecture

```
config.yaml ──▶ run.py ──▶ Snakemake 7 API ──▶ Snakefile
                                                   │ include
                                                   ▼
                                         workflow/rules/*.smk
                                                   │ shell, in each rule's conda env
                                                   ▼
                                         workflow/scripts/*
                                                   │
                                                   ▼
                          <output_dir>/results/<sample>/matchacell/…
```

`run.py` is a thin, matcha-themed wrapper around `snakemake.snakemake(**kwargs)`.
The `Snakefile` normalizes `output_dir`, includes one rule module per step and
per annotator, decides which optional annotators are enabled, and defines the
aggregating targets, which expand over every sample.

## Targets

| Target | What it builds |
| --- | --- |
| `cluster_stability` | Steps 0 and 1 for every sample: the clustered `.h5ad`, its `.rds` companion and the MatchA Verdict. |
| `annotation` | Step 2: the annotated `.h5ad` of every enabled annotator, running Step 1 first when needed. |

## The launcher: `run.py`

```
./run.py -w <target> -c config.yaml -q <cores> [options]
```

| Group | Flag | Purpose |
| --- | --- | --- |
| Required | `-w, --workflow` | Target(s) to build: `cluster_stability`, `annotation`. |
| | `-c, --configfile` | Path to the YAML config. |
| | `-q, --cores` | Cores available to Snakemake. |
| Execution | `-n, --dry-run` | Plan jobs without running them. |
| | `-f, --forceall` | Rebuild everything. |
| | `-ri, --rerun-incomplete` | Rerun jobs with incomplete outputs. |
| | `-u, --unlock` | Unlock the working directory. |
| | `--keep-going` | Continue independent jobs after a failure. |
| | `--latency-wait` | Seconds to wait for outputs on slow filesystems. |
| Reporting | `--list-workflows` | Show targets (including planned ones) and exit. |
| | `-d, --dag` | Emit the DAG in DOT format. |
| | `-p, --printshellcmds` | Print shell commands. |
| | `-l, --lint` | Run Snakemake lint (`text`/`json`). |
| Environment | `--no-conda` | Disable per-rule conda (on by default). |
| | `--conda-prefix` | Shared conda prefix directory. |
| | `--use-singularity` | Enable containers (off by default). |
| | `--bind`, `--singularity-args` | Container bind paths / extra args. |
| | `-di, --directory` | Working directory. |
| | `--snakefile` | Custom Snakefile path. |
| Advanced | `--resources` | `key=value` Snakemake resources, e.g. `mem_mb=64000`. |
| | `--allow-custom-target` | Run an arbitrary rule name. |

`run.py` warns if the active Snakemake is not 7.x. Singularity is **off by
default** because MatchACell ships conda environments rather than containers.
The annotator rules are per sample, so to run a single annotator, ask Snakemake
for its output file; see
[Running a single annotator](annotators/README.md#running-a-single-annotator).

## Configuration: `config.yaml`

```yaml
output_dir: "/path/to/output"     # results/ is created underneath
dtype: "h5ad"                     # "h5ad", or "rds" for Seurat/SingleCellExperiment input
type: "singlecell"

samples:
  pbmc3k: "data/pbmc3k.h5ad"        # sample_id: path-to-input

matchacell_cluster_stability:
  backend: "auto"                   # auto | gpu | cpu
  n_iter: 1000                      # bootstrap iterations per resolution
  extra: ""                         # extra flags forwarded to the script

matchacell_annotation:
  annot_file: "signatures.xlsx"     # marker-gene workbook (ScoreGenes, AddModuleScore, CIA)
  score_genes:    { ... }
  cia:            { ... }
  addmodulescore: { ... }
  celltypist:     { ... }
  scanvi:         { ... }           # skipped while reference_file is empty
  scparadise:     { ... }           # skipped while model_dir is empty
  cytetype:       { ... }           # skipped without an API token
  celltypeai:     { ... }           # skipped while tissue is empty
```

- **`samples`** maps a sample ID (used in output paths and report titles) to the
  path of its input. Add more `key: path` entries to fan out.
- **`dtype`** tells Step 0 whether inputs are `.h5ad` (copied as they are) or
  `.rds` (converted with `workflow/scripts/rds2h5.R`).
- **`backend`** chooses the Step 1 compute path. `auto` uses the GPU when
  `rapids-singlecell` imports successfully, otherwise CPU.
- **`n_iter`** is the bootstrap depth. More iterations smooth the Jaccard
  estimates at the cost of runtime.
- **`extra`** is forwarded verbatim to `matchacell_cluster_stability.py`, so any
  flag in [`cli-and-outputs.md`](cli-and-outputs.md) can be set there, e.g.
  `"--resolutions 0.2 0.5 1.0 2.0 --skip-tsne --use-hvg no"`. Use
  `--use-hvg no` before annotation; see
  [Keep every gene for annotation](annotators/README.md#keep-every-gene-for-annotation).
- **`matchacell_annotation`** holds one section per annotator, documented on each
  [annotator's page](annotators/README.md). Keep every section, even for
  annotators you ignore: most rules read theirs when the workflow is parsed.

## The rules

| Step | Rule | Module | Script | Environment |
| --- | --- | --- | --- | --- |
| 0 | `stage_h5` | `rds2h5.smk` | `rds2h5.R` for `.rds` input, otherwise a copy | `addmodulescore.yaml` |
| 1 | `matchacell_cluster_stability` | `cluster_stability.smk` | `matchacell_cluster_stability.py` | `matchacell.yaml` |
| 1 | `matchacell_cluster_stability_rds` | `cluster_stability.smk` | `h52rds.R` | `addmodulescore.yaml` |
| 2 | `score_genes` | `score_genes.smk` | `score_genes.py` | `matchacell.yaml` |
| 2 | `addmodulescore` | `addmodulescore.smk` | `addmodulescore_tool.R`, then `rds2h5.R` | `addmodulescore.yaml` |
| 2 | `cia` | `cia.smk` | `cia_tool.py` | `cia.yaml` |
| 2 | `celltypist` | `celltypist.smk` | `celltypist_tool.py` | `celltypist.yaml` |
| 2 | `scanvi` | `scanvi.smk` | `scanvi_tool.py` | `scanvi.yaml` |
| 2 | `scparadise` | `scparadise.smk` | `scparadise_tool.py` | `scparadise.yaml` |
| 2 | `cytetype` | `cytetype.smk` | `cytetype_tool.py` | `cytetype.yaml` |
| 2 | `celltypeai` | `celltypeai.smk` | `celltypeai_tool.py` | `celltypeai.yaml` |

Step 0 stages each sample as `results/<sample>/h5/<sample>.h5ad`, so later steps
never need to know the input format.

Step 1's rule reads that staged file and writes `clustered_multi_resolution.h5ad`
and `MatchA_Verdict.txt` (its tracked outputs) plus diagnostics under
`matchacell/`. A second rule converts the clustered object to a Seurat `.rds`
for R-side use, including AddModuleScore.

Each Step 2 rule reads the clustered object, takes the resolution from
`MatchA_Verdict.txt`, and writes its annotated `.h5ad` (its tracked output),
tables and plots to `annotation/<Annotator>/`. Because annotators are
independent rules, Snakemake runs them concurrently up to `--cores`.

## Adding an annotator

1. **Script**: `workflow/scripts/<tool>_tool.py` (or `.R`). Read the clustered
   `.h5ad`, get the resolution with `functions_annot.extract_best_res(verdict)`,
   and write `<tool>_annotated.h5ad` plus a per-cluster table to the output
   folder. `functions_annot.multi_umap` draws the per-label UMAP panels.
2. **Environment**: `workflow/envs/<tool>.yaml` with the tool's dependencies.
3. **Rule**: `workflow/rules/<tool>.smk`, following an existing annotator, and an
   `include:` line in the `Snakefile`.
4. **Config**: a section under `matchacell_annotation`. If the tool needs
   something users may not have (a model, a token, a server), read it with
   `.get()`, add a switch in the `Snakefile` (like `_SCANVI_REFERENCE`), and make
   its entry in `rule annotation` conditional on that switch, so a missing
   prerequisite skips the annotator instead of failing the run.
5. **Docs**: a page in `docs/annotators/`, a row in the
   [annotators overview](annotators/README.md), and a row in the README table.
6. **Test**: a dry-run scheduling test like `tests/test_scparadise_workflow.py`.

## Troubleshooting

- **`AttributeError: module 'snakemake' has no attribute 'snakemake'`**: you
  are on Snakemake 8.x. Install `snakemake-minimal=7.32.4` (see `environment.yml`).
- **`CreateCondaEnvironmentException … Non-conda folder exists at prefix`**: mamba
  2.x with Snakemake 7; see [installation](installation.md#mamba-2-and-snakemake-7).
- **Singularity errors**: leave Singularity off (the default); MatchACell uses
  conda environments.
- **`output_dir` produced a path like `.../outputresults/...`**: fixed in the
  `Snakefile` (the trailing slash is normalized), but double-check custom edits.
- **Conda environment paths**: each rule references `../envs/<name>.yaml`;
  repoint a rule's `conda:` directive if you keep a prebuilt environment
  elsewhere.
