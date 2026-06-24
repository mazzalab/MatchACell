"""Unit tests for argument parsing, theming and small helpers."""
from __future__ import annotations

import pytest

pytest.importorskip("scanpy")
import matchacell_cluster_stability as mc  # noqa: E402


def test_parse_args_defaults():
    args = mc.parse_args(["--input", "x.h5ad", "--outdir", "out"])
    assert args.backend == "auto"
    assert args.n_iter == 100
    assert args.fraction == pytest.approx(0.8)
    assert args.resolutions == [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]


def test_parse_args_overrides():
    args = mc.parse_args([
        "--input", "x.h5ad", "--outdir", "out",
        "--backend", "cpu", "--n-iter", "7",
        "--resolutions", "0.4", "0.9",
    ])
    assert args.backend == "cpu"
    assert args.n_iter == 7
    assert args.resolutions == [0.4, 0.9]


def test_sorted_leiden_keys_orders_numerically():
    ad = pytest.importorskip("anndata")
    import numpy as np

    a = ad.AnnData(np.zeros((3, 1), dtype="float32"))
    a.obs["leiden_1.0"] = ["0", "0", "0"]
    a.obs["leiden_0.2"] = ["0", "1", "0"]
    a.obs["leiden_0.05"] = ["0", "0", "1"]
    a.obs["not_leiden"] = ["a", "b", "c"]
    assert mc.sorted_leiden_keys(a) == ["leiden_0.05", "leiden_0.2", "leiden_1.0"]


def test_theme_respects_no_color(monkeypatch):
    mc.Theme.configure(no_color=True)
    assert mc.Theme.matcha("hello") == "hello"        # no ANSI when disabled


def test_matcha_scent_constant():
    assert "MatchACell" in mc.MATCHA_SCENT
