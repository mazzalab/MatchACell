# The pipeline

This document covers the Snakemake side: the launcher, the configuration, the
single rule that exists today, and how future annotator steps will slot in.

## Architecture

```
config.yaml ──▶ run.py ──▶ Snakemake 7 API ──▶ Snakefile
                                                   │ include
                                                   ▼
                                  workflow/rules/cluster_stability.smk
                                                   │ shell
                                                   ▼
                            workflow/scripts/matchacell_cluster_stability.py
                                                   │
                                                   ▼
                       <output_dir>/results/<sample>/matchacell/…
```

`run.py` is a thin, matcha-themed wrapper around `snakemake.snakemake(**kwargs)`.
The `Snakefile` normalizes `output_dir`, includes the Step-1 rule module, and
defines the aggregating target `cluster_stability` that expands over every
sample. The rule shells out to the Step-1 engine once per sample.

## The launcher: `run.py`

```
./run.py -w <target> -c config.yaml -q <cores> [options]
```

| Group | Flag | Purpose |
| --- | --- | --- |
| Required | `-w, --workflow` | Target(s) to build. Today: `cluster_stability`. |
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

## Configuration: `config.yaml`

```yaml
output_dir: "/path/to/output"     # results/ is created underneath
extension: "h5ad"

samples:
  pbmc3k: "data/pbmc3k.h5ad"        # sample_id: path-to-input-.h5ad

matchacell_cluster_stability:
  backend: "auto"                   # auto | gpu | cpu
  n_iter: 1000                      # bootstrap iterations per resolution
  extra: ""                         # extra flags forwarded to the script
```

- **`samples`** maps a sample ID (used in output paths and report titles) to the
  path of its input `.h5ad`. Add more `key: path` entries to fan out.
- **`backend`** chooses the compute path. `auto` uses the GPU when
  `rapids-singlecell` imports successfully, otherwise CPU.
- **`n_iter`** is the bootstrap depth. More iterations smooth the Jaccard
  estimates at the cost of runtime.
- **`extra`** is forwarded verbatim to `matchacell_cluster_stability.py`, so any
  flag in [`cli-and-outputs.md`](cli-and-outputs.md) can be set there, e.g.
  `"--resolutions 0.2 0.5 1.0 2.0 --skip-tsne --gpu-cell-threshold 50000"`.

## The rule

`workflow/rules/cluster_stability.smk` defines one rule:

```python
rule matchacell_cluster_stability:
    input:  h5ad   = lambda wc: config["samples"][wc.sample]
    output: clustered_h5ad = outputDir + "results/{sample}/matchacell/clustered_multi_resolution.h5ad"
    params: outdir, backend, n_iter, extra
    conda:  "../envs/matchacell.yaml"
    shell:  "python workflow/scripts/matchacell_cluster_stability.py --input … --outdir … --backend … --n-iter … {extra}"
```

The script writes `clustered_multi_resolution.h5ad` (the tracked output) plus a
tree of diagnostics under `params.outdir`.

## Adding future annotators (roadmap)

The repository is structured so additional steps drop in cleanly:

1. Add a script under `workflow/scripts/` (e.g. `annotate_<tool>.py`).
2. Add a rule module under `workflow/rules/` and `include:` it from the
   `Snakefile`.
3. Register a new target in `WORKFLOWS` in `run.py` (and remove it from
   `PLANNED`).
4. Run several annotators **in parallel** by giving Snakemake the cores and
   listing multiple targets: `./run.py -w annotate consensus -c config.yaml -q 32`.

Because annotators are independent rules, Snakemake schedules them concurrently
up to `--cores`.

## Troubleshooting

- **`AttributeError: module 'snakemake' has no attribute 'snakemake'`** — you
  are on Snakemake 8.x. Install `snakemake-minimal=7.32.4` (see `environment.yml`).
- **Singularity errors** — leave Singularity off (the default); MatchACell uses
  conda envs.
- **`output_dir` produced a path like `.../outputresults/...`** — fixed in the
  `Snakefile` (the trailing slash is normalized), but double-check custom edits.
- **Conda env path** — the rule references `../envs/matchacell.yaml`; repoint it
  if you keep a prebuilt environment elsewhere.
