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


def test_verdict_flags_largest_resolution_of_the_sweep():
    master = _master(
        [("0.2", str(c), 0.9) for c in range(2)]
        + [("0.5", str(c), 0.88) for c in range(4)]
    )
    v = mc.compute_verdict(master)
    assert v["recommended"] == "0.5" and v["at_sweep_edge"]
    assert "largest resolution swept" in mc.render_verdict_text(v)

    # a coarser resolution wins once the finest one is unstable
    master = _master(
        [("0.2", str(c), 0.9) for c in range(2)]
        + [("0.5", str(c), 0.5) for c in range(4)]
    )
    v = mc.compute_verdict(master)
    assert v["recommended"] == "0.2" and not v["at_sweep_edge"]


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


def test_verdict_counts_the_cells_in_unstable_clusters():
    # 0.5 passes the cluster-share bar (1 of 5 risky) but 30% of the cells sit
    # in that risky cluster: the coarser 0.2 wins. A tiny risky cluster passes.
    rows = ([("0.2", str(c), 0.9, 200) for c in range(3)]
            + [("0.5", str(c), 0.9, 100) for c in range(4)] + [("0.5", "4", 0.3, 180)])
    master = pd.DataFrame(rows, columns=["Resolution", "Cluster", "Jaccard_Stability", "n_cells"])
    v = mc.compute_verdict(master)
    assert v["recommended"] == "0.2"
    assert v["summary"].loc["0.5", "frac_risk_cells"] == pytest.approx(180 / 580)
    master.loc[master["Cluster"].eq("4"), "n_cells"] = 10
    v = mc.compute_verdict(master)
    assert v["recommended"] == "0.5" and v["frac_risk_cells"] == pytest.approx(10 / 410)
    assert "Cells in clusters < 0.60: 2.4%" in mc.render_verdict_text(v)


def test_summary_uses_the_given_thresholds():
    master = _master([("0.5", "0", 0.7), ("0.5", "1", 0.9)])
    s = mc.resolution_summary(master, high=0.95, risk=0.75)
    assert s.loc["0.5", "frac_high"] == 0 and s.loc["0.5", "frac_risk"] == 0.5
