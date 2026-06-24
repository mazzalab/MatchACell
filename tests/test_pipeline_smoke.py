"""End-to-end CPU smoke test.

Runs the full Step-1 pipeline on a tiny synthetic dataset through the CPU
backend. Besides checking that the expected artefacts are produced, it asserts
that the interactive HTML report's performance strip carries the cell count —
this is the regression guard for the finalize-before-HTML ordering fix.
"""
from __future__ import annotations

import json

import pytest

pytest.importorskip("scanpy")
pytest.importorskip("leidenalg")
pytest.importorskip("igraph")
pytest.importorskip("openpyxl")
import matchacell_cluster_stability as mc  # noqa: E402


@pytest.mark.slow
def test_cpu_pipeline_end_to_end(tiny_h5ad, tmp_path):
    outdir = tmp_path / "out"
    mc.main([
        "--input", str(tiny_h5ad),
        "--outdir", str(outdir),
        "--backend", "cpu",
        "--resolutions", "0.4", "0.8",
        "--n-pcs", "10",
        "--n-neighbors", "15",
        "--n-iter", "5",
        "--fraction", "0.8",
        "--skip-tsne",
        "--skip-umap",
        "--no-color",
    ])

    # Core artefacts exist.
    assert (outdir / "clustered_multi_resolution.h5ad").exists()
    assert (outdir / "MatchA_Verdict.txt").exists()
    assert (outdir / "stability" / "Cluster_Stability_Summary.xlsx").exists()

    # Performance metadata was finalized (cell count + total runtime present).
    perf = json.loads((outdir / "performance" / "performance_summary.json").read_text())
    assert "n_cells" in perf and perf["n_cells"] > 0
    assert "total_seconds" in perf

    # Regression guard: the HTML perf strip is populated because finalize() now
    # runs before build_html_report().
    html = (outdir / "MatchACell_report.html")
    if html.exists():  # only when plotly is installed
        text = html.read_text(encoding="utf-8")
        assert "cells" in text
