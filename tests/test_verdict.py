"""Unit tests for the MatchA Verdict logic."""
from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("scanpy")
import matchacell_cluster_stability as mc  # noqa: E402


def _master(rows):
    return pd.DataFrame(rows, columns=["Resolution", "Cluster", "Jaccard_Stability"])


def test_resolution_summary_columns():
    master = _master([
        ("0.2", "0", 0.9), ("0.2", "1", 0.8),
        ("0.5", "0", 0.7), ("0.5", "1", 0.5),
    ])
    s = mc.resolution_summary(master)
    assert list(s.index) == ["0.2", "0.5"]            # sorted numerically
    for col in ("n_clusters", "median", "mean", "min", "frac_high", "frac_risk"):
        assert col in s.columns
    assert s.loc["0.2", "n_clusters"] == 2


def test_verdict_picks_finest_stable_resolution():
    # 0.2: 2 stable clusters; 0.5: 4 stable clusters -> prefer the finer 0.5.
    master = _master(
        [("0.2", str(c), 0.9) for c in range(2)]
        + [("0.5", str(c), 0.88) for c in range(4)]
    )
    v = mc.compute_verdict(master)
    assert v["recommended"] == "0.5"
    assert v["n_clusters"] == 4
    assert "finest" in v["basis"]
    assert len(v["unstable_clusters"]) == 0


def test_verdict_falls_back_when_nothing_qualifies():
    # Every resolution is risky -> fallback to the highest median.
    master = _master([
        ("0.2", "0", 0.50), ("0.2", "1", 0.40),
        ("0.5", "0", 0.30), ("0.5", "1", 0.20),
    ])
    v = mc.compute_verdict(master)
    assert v["recommended"] == "0.2"                  # higher median
    assert "fallback" in v["basis"]


def test_render_verdict_text_contains_scent_and_recommendation():
    master = _master([("0.3", "0", 0.9), ("0.3", "1", 0.86)])
    v = mc.compute_verdict(master)
    text = mc.render_verdict_text(v)
    assert "MatchA Verdict" in text
    assert "leiden_0.3" in text
    assert mc.MATCHA_SCENT in text
