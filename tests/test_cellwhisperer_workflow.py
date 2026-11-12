"""Validate CellWhisperer scheduling with the real Snakemake parser, without running models."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake")
pytestmark = pytest.mark.skipif(SNAKEMAKE is None, reason="Snakemake is not installed")


@pytest.fixture
def workflow_config(tmp_path):
    config = yaml.safe_load((ROOT / "config.yaml").read_text())
    config["output_dir"] = str(tmp_path / "output with spaces")
    sample = tmp_path / "input.h5ad"
    sample.touch()
    config["samples"] = {"sample1": str(sample)}
    signatures = tmp_path / "signatures.xlsx"
    signatures.touch()
    config["matchacell_annotation"]["annot_file"] = str(signatures)
    config["matchacell_annotation"]["celltypeai"]["tissue"] = ""
    # Satisfy the other annotators to isolate the new branch of the DAG.
    base = Path(config["output_dir"]) / "results/sample1/matchacell"
    for method, filename in [
        ("ScoreGenes", "score_genes"), ("CIA", "cia"),
        ("CellTypist", "celltypist"), ("AddModuleScore", "addmodulescore"),
    ]:
        output = base / "annotation" / method / f"{filename}_annotated.h5ad"
        output.parent.mkdir(parents=True)
        output.touch()
    return config


@pytest.fixture
def cw_repo(tmp_path):
    """A fake CellWhisperer checkout: a pixi manifest and the default checkpoint."""
    repo = tmp_path / "cell whisperer"
    ckpt = repo / "results/models/jointemb/cellwhisperer_clip_v1.ckpt"
    ckpt.parent.mkdir(parents=True)
    ckpt.touch()
    (repo / "pixi.toml").touch()
    return repo


def dry_run(config, tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    env = os.environ.copy()
    env.pop("CYTETYPE_API_TOKEN", None)
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    return subprocess.run(
        [SNAKEMAKE, "--snakefile", str(ROOT / "Snakefile"), "--configfile", str(path),
         "--cores", "2", "--dry-run", "--printshellcmds", "annotation"],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=60,
    )


@pytest.mark.parametrize("omit_section", [False, True])
def test_unconfigured_cellwhisperer_is_skipped(workflow_config, tmp_path, omit_section):
    if omit_section:
        del workflow_config["matchacell_annotation"]["cellwhisperer"]
    result = dry_run(workflow_config, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r"^cellwhisperer\s+\d+", result.stdout, re.MULTILINE) is None


def test_configured_cellwhisperer_schedules_with_workbook_labels(workflow_config, tmp_path, cw_repo):
    workflow_config["matchacell_annotation"]["cellwhisperer"].update(
        repo_dir=str(cw_repo), template="{label} from human pancreas",
    )
    result = dry_run(workflow_config, tmp_path)
    out = result.stdout
    assert result.returncode == 0, out + result.stderr
    assert re.search(r"^cellwhisperer\s+1", out, re.MULTILINE)
    assert re.search(r"^matchacell_cluster_stability\s+1", out, re.MULTILINE)
    assert f"--manifest-path '{cw_repo / 'pixi.toml'}'" in out
    assert "cellwhisperer_clip_v1.ckpt" in out
    assert "--template '{label} from human pancreas'" in out
    assert "--n_threads 2" in out
    assert "--labels" not in out  # candidates come from the workbook sheets


def test_explicit_labels_are_passed(workflow_config, tmp_path, cw_repo):
    workflow_config["matchacell_annotation"]["cellwhisperer"].update(
        repo_dir=str(cw_repo), labels=["beta cell", "TFF1+ ductal cell"],
    )
    result = dry_run(workflow_config, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "--labels 'beta cell' 'TFF1+ ductal cell'" in result.stdout


def test_missing_checkpoint_fails(workflow_config, tmp_path, cw_repo):
    (cw_repo / "results/models/jointemb/cellwhisperer_clip_v1.ckpt").unlink()
    workflow_config["matchacell_annotation"]["cellwhisperer"]["repo_dir"] = str(cw_repo)
    result = dry_run(workflow_config, tmp_path)
    assert result.returncode != 0
    assert "cellwhisperer_clip_v1.ckpt" in result.stdout + result.stderr
