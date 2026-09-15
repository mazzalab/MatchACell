"""run.py: PBS submission, per-config locking and per-rule cluster resources."""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

pytest.importorskip("snakemake")

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("matchacell_run", ROOT / "run.py")
run = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run)


def test_cluster_cmd_uses_pbs_select_syntax(tmp_path):
    cmd = run.build_cluster_cmd("workq", tmp_path)
    assert cmd.startswith("qsub -q workq ")
    assert "-l select=1:ncpus={threads}:mem={resources.mem_mb}mb" in cmd
    assert "-l walltime={resources.runtime}:00" in cmd
    assert "-N smk.{rule}.{jobid}" in cmd
    assert f"-o {tmp_path}/{{rule}}.{{jobid}}.out" in cmd
    assert f"-e {tmp_path}/{{rule}}.{{jobid}}.err" in cmd


def test_cluster_defaults():
    args = run.build_parser().parse_args(["-w", "annotation", "-c", "x.yaml", "-cl"])
    assert (args.cluster, args.queue, args.jobs, args.restart_times) == (True, "workq", 20, 2)
    assert args.conda_frontend == "conda"
    local = run.build_parser().parse_args(["-w", "annotation", "-c", "x.yaml"])
    assert not local.cluster


def test_same_config_cannot_run_twice(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "THIS_DIR", tmp_path)
    config, other = tmp_path / "config.yaml", tmp_path / "other.yaml"
    config.touch()
    other.touch()
    first = run.acquire_project_lock(config)
    with pytest.raises(SystemExit):
        run.acquire_project_lock(config)
    second = run.acquire_project_lock(other)   # a different project runs alongside
    first.close()                               # the lock goes with the process/handle
    run.acquire_project_lock(config).close()
    second.close()


def _rule_blocks():
    for module in sorted((ROOT / "workflow" / "rules").glob("*.smk")):
        text = module.read_text()
        for match in re.finditer(r"^rule (\w+):\n(.*?)(?=^rule |\Z)", text, re.MULTILINE | re.DOTALL):
            yield module.name, match.group(1), match.group(2)


def test_every_submitted_rule_requests_memory_and_walltime():
    # Snakemake prints a rule's `message:` instead of its resources, so they
    # can't be read back from a dry-run; the real check is `qstat -xf` on a run.
    blocks = list(_rule_blocks())
    assert len(blocks) >= 11
    for module, name, body in blocks:
        for key in ("mem_mb", "runtime"):
            assert re.search(rf'{key}=cluster_resource\("{name}", "{key}", \d+\)', body), (
                f"{module}: rule {name} has no {key} for PBS submission")
