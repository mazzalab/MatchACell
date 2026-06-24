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
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

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
    print(f"{matcha('Cores:')}            {args.cores}")
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
    print(f"{matcha('Use conda:')}        {not args.no_conda}")
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
    print()


def validate_args(args, snakefile: Path, configfile: Path) -> None:
    if not snakefile.exists():
        die(f"Cannot find Snakefile: {snakefile}")
    if not configfile.exists():
        die(f"Cannot find config file: {configfile}")
    if not args.unlock and not args.workflow:
        die("No workflow target provided. Use -w cluster_stability (see --list-workflows).")
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

    print_run_summary(args, snakefile, configfile, config, singularity_args, resources)

    snakemake_kwargs: Dict[str, Any] = {
        "snakefile": str(snakefile),
        "configfiles": [str(configfile)],
        "targets": args.workflow or [],
        "workdir": workdir,
        "cores": args.cores,
        "dryrun": args.dry_run,
        "use_conda": not args.no_conda,
        "use_singularity": args.use_singularity,
        "singularity_args": singularity_args,
        "forceall": args.forceall,
        "force_incomplete": args.rerun_incomplete,
        "unlock": args.unlock,
        "printdag": args.dag,
        "lint": args.lint,
        "printshellcmds": args.printshellcmds,
        "keepgoing": args.keep_going,
        "latency_wait": args.latency_wait,
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

  {matcha('./run.py -w cluster_stability -c config.yaml -q 16')}
      Run Step 1.

  {matcha('./run.py --list-workflows')}
      Show available workflow targets.

  {matcha('./run.py -w cluster_stability -c config.yaml -q 16 --dag | dot -Tsvg > dag.svg')}
      Export the DAG as SVG.
""",
    )

    required = parser.add_argument_group("Required for normal execution")
    required.add_argument("-w", "--workflow", nargs="+",
                          help="Workflow target(s). Currently: cluster_stability.")
    required.add_argument("-c", "--configfile",
                          help="Path to the YAML/JSON Snakemake config file.")
    required.add_argument("-q", "--cores", type=int, default=1,
                          help="Number of CPU cores available to Snakemake.")

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
                          help="Allow running rule names not listed among the main targets.")

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
