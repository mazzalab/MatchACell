#!/usr/bin/env python3
"""
MatchACell — matcha-themed launcher for the Snakemake workflow.

This wraps the Snakemake **7.x** Python API (``snakemake.snakemake(**kwargs)``).
It is intentionally pinned to Snakemake 7 (developed against 7.32.4); the 8.x
API is incompatible.

Examples
--------
Dry-run Step 1 (QC + cluster-stability):
    ./run.py -w cluster_stability -c config.yaml -q 16 -n

Run Step 1:
    ./run.py -w cluster_stability -c config.yaml -q 16

Run Step 2 (every enabled annotator):
    ./run.py -w annotation -c config.yaml -q 16

Submit each job as its own PBS job, at most 20 at a time:
    ./run.py -w annotation -c config.yaml -q 1 -cl -qu workq -j 20

Build the per-rule conda environments without running anything:
    ./run.py -w annotation -c config.yaml -q 4 --conda-create-envs-only

Run with shell commands and reasons printed:
    ./run.py -w cluster_stability -c config.yaml -q 16 --printshellcmds

List available workflow targets:
    ./run.py --list-workflows

Export the DAG as SVG:
    ./run.py -w cluster_stability -c config.yaml -q 16 --dag | dot -Tsvg > dag.svg

Unlock the working directory after an interrupted run:
    ./run.py -c config.yaml -q 1 --unlock
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import os
import platform
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# On Apple Silicon, conda's solver can silently resolve some per-rule envs
# onto osx-64 (Rosetta-emulated) builds instead of native osx-arm64 -- e.g.
# scvi.yaml's jax/jaxlib pulled an x86_64 wheel this way and crashed with
# "This version of jaxlib was built using AVX instructions, which your CPU
# ... do not support" (2026-09-07, workflow/envs/scanvi.yaml). Pin the
# subdir so `conda env create` is constrained to arm64-native packages (and
# fails loudly if one genuinely isn't available, instead of quietly falling
# back to an emulated build). Only set if the caller hasn't already chosen
# one explicitly.
if sys.platform == "darwin" and platform.machine() == "arm64":
    os.environ.setdefault("CONDA_SUBDIR", "osx-arm64")

import snakemake

try:  # PyYAML only needed for the run summary
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

try:  # colorama only helps legacy Windows terminals; truecolor is used below
    from colorama import just_fix_windows_console

    just_fix_windows_console()
except Exception:  # pragma: no cover
    pass


THIS_DIR = Path(__file__).resolve().parent

# PBS fallback for rules that don't declare `resources: mem_mb`/`runtime`
# (every MatchACell rule does; see cluster_resource() in the Snakefile).
# Only used with --cluster; local runs are unaffected.
CLUSTER_DEFAULT_MEM_MB = 8000
CLUSTER_DEFAULT_RUNTIME_MIN = 120

# --------------------------------------------------------------------------- #
# Matcha theme (truecolor; mirrors workflow/scripts/matchacell_cluster_stability.py)
# --------------------------------------------------------------------------- #
_MATCHA = (122, 182, 97)   # bright matcha leaf
_DEEP = (74, 124, 60)      # steeped deep green
_LIGHT = (183, 213, 160)   # foam green
_STONE = (140, 150, 120)   # muted sage
_WARN = (200, 170, 60)     # amber
_ERR = (200, 90, 60)       # burnt

MATCHA_SCENT = "Brewed with MatchACell — matcha-grade single-cell consensus"


def _color_enabled() -> bool:
    return os.environ.get("NO_COLOR") is None and sys.stdout.isatty()


def paint(text: str, rgb=_MATCHA, bold: bool = False) -> str:
    if not _color_enabled():
        return text
    b = "\033[1m" if bold else ""
    return f"{b}\033[38;2;{rgb[0]};{rgb[1]};{rgb[2]}m{text}\033[0m"


def matcha(t, bold=False):
    return paint(t, _MATCHA, bold)


def deep(t, bold=False):
    return paint(t, _DEEP, bold)


def dim(t):
    return paint(t, _STONE)


def die(message: str, exit_code: int = 1) -> None:
    sys.stderr.write(paint(f"ERROR: {message}\n", _ERR, bold=True))
    raise SystemExit(exit_code)


def warn(message: str) -> None:
    sys.stderr.write(paint(f"WARNING: {message}\n", _WARN, bold=True))


def info(message: str) -> None:
    sys.stderr.write(matcha(f"{message}\n"))


def success(message: str) -> None:
    sys.stderr.write(paint(f"{message}\n", _DEEP, bold=True))


# --------------------------------------------------------------------------- #
# Branding
# --------------------------------------------------------------------------- #
LOGO = r"""
        )  )  (
       (  (  ) )
        )  ) ( (
      .-~~~~~~~-.
     ( MatchACell )
      `._______.'~
       |       |    {tagline}
       |       |
        \_____/
"""


WORKFLOWS: Dict[str, str] = {
    "cluster_stability": (
        "Step 1 — data-driven QC + multi-resolution Leiden cluster-stability "
        "optimizer, producing the MatchA Verdict."
    ),
    "annotation": (
        "Step 2 — cell type annotation with every enabled annotator."
    ),
}

# Targets that are on the roadmap but not yet wired as rules. Shown (dimmed) by
# --list-workflows so the multi-annotator direction is visible.
PLANNED: Dict[str, str] = {
    "annotate": "Step 2 — multiple cell type/state annotators run in parallel (planned).",
    "consensus": "Step 3 — cross-annotator consensus + unified report (planned).",
}


def print_banner() -> None:
    print(matcha(LOGO.format(tagline=deep("single-cell consensus")), bold=True))
    print(deep("  Multi-Annotator Consensus for single-cell types & states"))
    print(dim("  " + MATCHA_SCENT))
    print()


def print_workflows() -> None:
    print_banner()
    print(matcha("Available workflow targets:", bold=True))
    print()
    for name, description in WORKFLOWS.items():
        print(f"  {matcha(name, bold=True):<24} {description}")
    print()
    print(dim("Planned (not yet wired):"))
    for name, description in PLANNED.items():
        print(dim(f"  {name:<22} {description}"))
    print()


# --------------------------------------------------------------------------- #
# Config helpers
# --------------------------------------------------------------------------- #
def load_config(configfile: Path) -> Dict[str, Any]:
    if yaml is None:
        warn("PyYAML is not installed; config summary will be skipped.")
        return {}
    try:
        with configfile.open("r") as handle:
            data = yaml.safe_load(handle)
    except Exception as exc:
        warn(f"Could not parse config file for summary: {exc}")
        return {}
    return data or {}


def normalize_samples(samples_obj: Any) -> List[str]:
    if samples_obj is None:
        return []
    if isinstance(samples_obj, dict):
        return list(samples_obj.keys())
    if isinstance(samples_obj, list):
        return [str(x) for x in samples_obj]
    return [str(samples_obj)]


def parse_resources(resource_args: Optional[List[str]]) -> Optional[Dict[str, int]]:
    if not resource_args:
        return None
    resources: Dict[str, int] = {}
    for item in resource_args:
        if "=" not in item:
            die(f"Invalid resource format: {item}. Use key=value, e.g. mem_mb=64000")
        key, value = item.split("=", 1)
        key, value = key.strip(), value.strip()
        if not key:
            die(f"Invalid resource key in: {item}")
        try:
            resources[key] = int(value)
        except ValueError:
            die(f"Resource value must be an integer in: {item}")
    return resources


def build_singularity_args(args: argparse.Namespace) -> str:
    pieces: List[str] = []
    if args.singularity_args:
        pieces.append(args.singularity_args.strip())
    for bind_path in args.bind:
        pieces.append(f"--bind {bind_path}")
    return " ".join(pieces).strip()


# --------------------------------------------------------------------------- #
# Locking and PBS submission
# --------------------------------------------------------------------------- #
def acquire_project_lock(configfile: Path):
    """Refuse a second concurrent run of the same config, without locking the repo.

    Snakemake's own lock covers the whole working directory, so runs of two
    *different* configs from this repository would wait on each other. This lock
    is keyed on the config file's absolute path instead (Snakemake's is disabled
    with lock=False), and the kernel releases the flock however the process exits,
    so a killed run needs no --unlock. The caller must keep the returned handle
    referenced for the life of the process.
    """
    lock_dir = THIS_DIR / ".smk_locks"
    lock_dir.mkdir(exist_ok=True)
    key = hashlib.sha1(str(configfile.resolve()).encode()).hexdigest()[:16]
    lock_path = lock_dir / f"{key}.lock"
    handle = open(lock_path, "a+")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        handle.close()
        die(
            f"Another run.py is already running config '{configfile}' (lock: {lock_path}). "
            "Two runs of the same config would race on the same output files; wait for "
            "the other one to finish."
        )
    handle.seek(0)
    handle.truncate()
    handle.write(str(os.getpid()))
    handle.flush()
    return handle


def build_cluster_cmd(queue: str, logdir: Path) -> str:
    """qsub command Snakemake runs for each job (OpenPBS / PBS Pro `select` syntax).

    Doubled braces are Snakemake's own per-job placeholders and survive this
    .format() call; only the queue and log directory are filled in here.
    """
    template = (
        "qsub -q {queue} "
        "-l select=1:ncpus={{threads}}:mem={{resources.mem_mb}}mb "
        "-l walltime={{resources.runtime}}:00 "
        "-N smk.{{rule}}.{{jobid}} "
        "-o {logdir}/{{rule}}.{{jobid}}.out "
        "-e {logdir}/{{rule}}.{{jobid}}.err"
    )
    return template.format(queue=queue, logdir=logdir)


def print_run_summary(args, snakefile, configfile, config, singularity_args, resources) -> None:
    samples = normalize_samples(config.get("samples"))
    output_dir = config.get("output_dir", "not found in config")
    mc = config.get("matchacell_cluster_stability", {}) or {}

    print_banner()
    print(matcha("Run summary", bold=True))
    print(matcha("-----------"))
    print(f"{matcha('Snakefile:')}        {snakefile}")
    print(f"{matcha('Config:')}           {configfile}")
    print(f"{matcha('Targets:')}          {', '.join(args.workflow or [])}")
    cores_note = " (ignored with --cluster)" if args.cluster else ""
    print(f"{matcha('Cores:')}            {args.cores}{cores_note}")
    print(f"{matcha('Workdir:')}          {args.directory or os.getcwd()}")
    print(f"{matcha('Output dir:')}       {output_dir}")
    print(f"{matcha('Backend:')}          {mc.get('backend', 'n/a')}")
    print(f"{matcha('Bootstrap iters:')}  {mc.get('n_iter', 'n/a')}")
    print(f"{matcha('Samples:')}          {len(samples)}")
    if samples:
        preview = ", ".join(samples[:8])
        if len(samples) > 8:
            preview += f", … +{len(samples) - 8} more"
        print(f"{matcha('Sample IDs:')}       {preview}")
    print(f"{matcha('Dry run:')}          {args.dry_run}")
    print(f"{matcha('Use conda:')}        {not args.no_conda} (frontend: {args.conda_frontend})")
    if args.conda_create_envs_only:
        print(matcha("Conda envs only:   yes (nothing else runs)"))
    print(f"{matcha('Use Singularity:')}  {args.use_singularity}")
    if args.use_singularity and singularity_args:
        print(f"{matcha('Singularity args:')} {singularity_args}")
    if resources:
        print(f"{matcha('Resources:')}        {resources}")
    if args.forceall:
        print(matcha("Force all jobs:    yes"))
    if args.rerun_incomplete:
        print(matcha("Rerun incomplete:  yes"))
    if args.keep_going:
        print(matcha("Keep going:        yes"))
    if args.cluster:
        print(f"{matcha('PBS cluster:')}      queue {args.queue}, at most {args.jobs} jobs at once")
        print(f"{matcha('Restart times:')}    {args.restart_times}")
        print(f"{matcha('PBS logs:')}         {THIS_DIR / 'logs' / 'pbs'}")
    print()


def validate_args(args, snakefile: Path, configfile: Path) -> None:
    if not snakefile.exists():
        die(f"Cannot find Snakefile: {snakefile}")
    if not configfile.exists():
        die(f"Cannot find config file: {configfile}")
    if not args.unlock and not args.workflow:
        die(f"No workflow target provided. Use -w {' or '.join(WORKFLOWS)} (see --list-workflows).")
    if args.workflow and not args.allow_custom_target:
        known = {**WORKFLOWS, **PLANNED}
        invalid = [t for t in args.workflow if t not in known]
        if invalid:
            valid = ", ".join(WORKFLOWS.keys())
            die(
                f"Unknown workflow target(s): {', '.join(invalid)}. "
                f"Available now: {valid}. "
                "Use --allow-custom-target to run any rule name manually."
            )
        planned = [t for t in args.workflow if t in PLANNED]
        if planned:
            die(f"Target(s) {', '.join(planned)} are planned but not yet wired as rules.")


def _check_snakemake_version() -> None:
    try:
        major = int(snakemake.__version__.split(".")[0])
    except Exception:  # pragma: no cover
        return
    if major != 7:
        warn(
            f"run.py targets the Snakemake 7.x Python API but detected "
            f"{snakemake.__version__}. Pin snakemake-minimal=7.32.4 "
            "(see environment.yml)."
        )


def run_snakemake(args: argparse.Namespace) -> int:
    _check_snakemake_version()
    snakefile = Path(args.snakefile).resolve() if args.snakefile else THIS_DIR / "Snakefile"
    configfile = Path(args.configfile).resolve()
    workdir = str(Path(args.directory).resolve()) if args.directory else None

    validate_args(args, snakefile, configfile)

    config = load_config(configfile)
    resources = parse_resources(args.resources)
    singularity_args = build_singularity_args(args)

    # Held for the whole run; replaces Snakemake's directory-wide lock (lock=False).
    _project_lock = acquire_project_lock(configfile)  # noqa: F841

    cores = args.cores
    cluster_kwargs: Dict[str, Any] = {}
    if args.cluster:
        logdir = THIS_DIR / "logs" / "pbs"
        logdir.mkdir(parents=True, exist_ok=True)
        from snakemake.resources import DefaultResources

        cluster_kwargs = {
            "cluster": build_cluster_cmd(args.queue, logdir),
            "nodes": args.jobs,
            "restart_times": args.restart_times,
            "default_resources": DefaultResources(
                [f"mem_mb={CLUSTER_DEFAULT_MEM_MB}", f"runtime={CLUSTER_DEFAULT_RUNTIME_MIN}"]
            ),
        }
        # Snakemake caps each rule's threads at `cores` before filling {threads}
        # into the qsub template, even in cluster mode: with -q 1, an 8-thread rule
        # would be submitted with ncpus=1. A sentinel above any rule's threads keeps
        # them intact; -q itself doesn't matter once --cluster is set.
        cores = 999

    print_run_summary(args, snakefile, configfile, config, singularity_args, resources)

    snakemake_kwargs: Dict[str, Any] = {
        "snakefile": str(snakefile),
        "configfiles": [str(configfile)],
        "targets": args.workflow or [],
        "workdir": workdir,
        "cores": cores,
        "dryrun": args.dry_run,
        "use_conda": not args.no_conda,
        "conda_frontend": args.conda_frontend,
        "conda_create_envs_only": args.conda_create_envs_only,
        "use_singularity": args.use_singularity,
        "singularity_args": singularity_args,
        "forceall": args.forceall,
        "force_incomplete": args.rerun_incomplete,
        "unlock": args.unlock,
        "lock": False,
        "printdag": args.dag,
        "lint": args.lint,
        "printshellcmds": args.printshellcmds,
        "keepgoing": args.keep_going,
        "latency_wait": args.latency_wait,
        **cluster_kwargs,
    }
    if resources:
        snakemake_kwargs["resources"] = resources
    if args.wrapper_prefix:
        snakemake_kwargs["wrapper_prefix"] = args.wrapper_prefix
    if args.conda_prefix:
        snakemake_kwargs["conda_prefix"] = args.conda_prefix

    status = snakemake.snakemake(**snakemake_kwargs)

    if status:
        success("Workflow finished successfully. 🍵".replace(" 🍵", ""))
        return 0
    die("Workflow failed.", exit_code=1)
    return 1  # unreachable; keeps type checkers happy


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=matcha("Matcha-themed launcher for the MatchACell Snakemake workflow."),
        epilog=f"""
Examples:

  {matcha('./run.py -w cluster_stability -c config.yaml -q 16 -n')}
      Dry-run Step 1 (QC + cluster stability).

  {matcha('./run.py -w annotation -c config.yaml -q 16')}
      Run Step 2 (and Step 1 first, if needed) on this machine.

  {matcha('./run.py -w annotation -c config.yaml -q 1 -cl -qu workq -j 20')}
      Submit each job as its own PBS job via qsub, at most 20 at a time.

  {matcha('./run.py -w annotation -c config.yaml -q 4 --conda-create-envs-only')}
      Only build the per-rule conda environments.

  {matcha('./run.py --list-workflows')}
      Show available workflow targets.

  {matcha('./run.py -w cluster_stability -c config.yaml -q 16 --dag | dot -Tsvg > dag.svg')}
      Export the DAG as SVG.
""",
    )

    required = parser.add_argument_group("Required for normal execution")
    required.add_argument("-w", "--workflow", nargs="+",
                          help=f"Workflow target(s): {', '.join(WORKFLOWS)}.")
    required.add_argument("-c", "--configfile",
                          help="Path to the YAML/JSON Snakemake config file.")
    required.add_argument("-q", "--cores", type=int, default=1,
                          help="Number of CPU cores available to Snakemake (ignored with --cluster).")

    execution = parser.add_argument_group("Execution mode")
    execution.add_argument("-n", "--dry-run", action="store_true",
                           help="Perform a dry-run without executing jobs.")
    execution.add_argument("-f", "--forceall", action="store_true",
                           help="Force execution of all jobs even if outputs exist.")
    execution.add_argument("-ri", "--rerun-incomplete", action="store_true",
                           help="Rerun jobs with incomplete outputs.")
    execution.add_argument("-u", "--unlock", action="store_true",
                           help="Unlock the working directory after an interrupted run.")
    execution.add_argument("--keep-going", action="store_true",
                           help="Continue independent jobs after one job fails.")
    execution.add_argument("--latency-wait", type=int, default=60,
                           help="Seconds to wait for outputs on slow filesystems.")

    cluster = parser.add_argument_group("Cluster (PBS)")
    cluster.add_argument("-cl", "--cluster", action="store_true",
                         help="Submit each rule instance as its own PBS job via qsub, instead of "
                              "running everything locally (default off; plain -q N still runs "
                              "everything locally with N cores).")
    cluster.add_argument("-qu", "--queue", type=str, default="workq",
                         help="PBS queue to submit to when --cluster is set.")
    cluster.add_argument("-j", "--jobs", type=int, default=20,
                         help="Max number of concurrent PBS jobs when --cluster is set.")
    cluster.add_argument("-rt", "--restart-times", type=int, default=2,
                         help="With --cluster, resubmit a failed job up to this many times; each "
                              "attempt multiplies the rule's mem_mb and runtime by the attempt number.")

    reporting = parser.add_argument_group("Reporting and debugging")
    reporting.add_argument("--list-workflows", "--show-workflows",
                           dest="show_workflows", action="store_true",
                           help="Print available workflow targets and exit.")
    reporting.add_argument("-d", "--dag", action="store_true",
                           help="Print DAG in DOT format (pipe to `dot -Tsvg`).")
    reporting.add_argument("--printshellcmds", "-p", action="store_true",
                           help="Print shell commands before executing them.")
    reporting.add_argument("-l", "--lint", default=None,
                           help='Run Snakemake lint. Common values: "text" or "json".')

    environment = parser.add_argument_group("Environment and paths")
    environment.add_argument("-di", "--directory", type=str, default=None,
                             help="Working directory for Snakemake execution.")
    environment.add_argument("--snakefile", type=str, default=None,
                             help="Custom Snakefile path. Default: Snakefile next to run.py.")
    environment.add_argument("--no-conda", action="store_true",
                             help="Disable Snakemake conda integration (on by default).")
    environment.add_argument("--conda-prefix", type=str, default=None,
                             help="Optional Snakemake conda prefix directory.")
    environment.add_argument("--conda-frontend", choices=["conda", "mamba"], default="conda",
                             help="Tool that creates the per-rule environments. Snakemake 7 "
                                  "fails to create them with mamba 2.x.")
    environment.add_argument("--conda-create-envs-only", action="store_true",
                             help="Only create the conda environments the targets need, then exit.")
    environment.add_argument("--use-singularity", action="store_true",
                             help="Enable Singularity/Apptainer (off by default; "
                                  "MatchACell uses conda envs).")
    environment.add_argument("--bind", action="append", default=["/data"],
                             help="Singularity bind path. Can be used multiple times.")
    environment.add_argument("--singularity-args", type=str, default="",
                             help="Additional raw Singularity arguments.")
    environment.add_argument("-wr", "--wrapper-prefix", type=str, default=None,
                             help="Prefix for a local Snakemake wrapper repository.")

    advanced = parser.add_argument_group("Advanced")
    advanced.add_argument("--resources", nargs="+", default=None,
                          help="Custom Snakemake resources, e.g. --resources mem_mb=64000.")
    advanced.add_argument("--allow-custom-target", action="store_true",
                          help="Allow running rule names or output files not listed among the "
                               "main targets.")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.show_workflows:
        print_workflows()
        return 0

    if not args.configfile:
        die("Missing config file. Use -c config.yaml")

    return run_snakemake(args)


if __name__ == "__main__":
    raise SystemExit(main())
