"""Validate scheduling with the real Snakemake parser, without running models."""
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
def test_unconfigured_scparadise_is_skipped(workflow_config, tmp_path, omit_section):
    if omit_section:
        del workflow_config["matchacell_annotation"]["scparadise"]
    result = dry_run(workflow_config, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r"^scparadise\s+\d+", result.stdout, re.MULTILINE) is None


@pytest.mark.parametrize("detect_unknown", [False, True])
def test_configured_scparadise_schedules_dependencies(workflow_config, tmp_path, detect_unknown):
    model = tmp_path / "models with spaces" / "Human_PBMC_scAdam"
    model.mkdir(parents=True)
    (model / "model_v2.pth").touch()
    if detect_unknown:
        (model / "unknown_detector.json").touch()
    workflow_config["matchacell_annotation"]["scparadise"].update(
        model_dir=str(model), detect_unknown=detect_unknown, layer="",
    )
    result = dry_run(workflow_config, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r"^scparadise\s+1", result.stdout, re.MULTILINE)
    assert re.search(r"^matchacell_cluster_stability\s+1", result.stdout, re.MULTILINE)
    assert "MatchA_Verdict.txt" in result.stdout
    assert f"--model_dir '{model}'" in result.stdout
    assert "--n_threads 2" in result.stdout
    assert "--layer=" in result.stdout
    assert ("--detect_unknown" in result.stdout) == detect_unknown


def test_missing_configured_model_fails(workflow_config, tmp_path):
    workflow_config["matchacell_annotation"]["scparadise"]["model_dir"] = str(tmp_path / "missing_scAdam")
    result = dry_run(workflow_config, tmp_path)
    assert result.returncode != 0
    assert "model_v2.pth" in result.stdout + result.stderr
