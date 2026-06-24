"""Unit tests for the data-driven QC helpers."""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scanpy")
import matchacell_cluster_stability as mc  # noqa: E402


def test_is_outlier_flags_both_tails():
    # A spread bulk (nonzero MAD) plus two extreme values in either tail.
    vals = np.concatenate([np.arange(5, 16, dtype=float), [1000.0, -500.0]])
    mask = mc.is_outlier(vals, nmads=3.0)
    assert mask[-2] and mask[-1]          # 1000 and -500 are outliers
    assert not mask[:-2].any()            # the bulk is not


def test_is_outlier_zero_mad_returns_all_false():
    vals = np.full(20, 7.0)
    assert not mc.is_outlier(vals, nmads=3.0).any()


def test_mad_bounds_symmetry():
    vals = np.arange(100, dtype=float)
    lo, hi = mc.mad_bounds(vals, nmads=2.0)
    med = np.median(vals)
    assert lo < med < hi
    assert pytest.approx(hi - med, rel=1e-9) == med - lo


def test_flag_feature_classes(tiny_adata):
    flags = mc.flag_feature_classes(tiny_adata)
    assert flags["mt"] == 2          # MT-CO1, MT-ND1
    assert flags["control"] == 1     # NegControlProbe_*
    assert tiny_adata.var["mt"].sum() == 2
    assert tiny_adata.var["control"].sum() == 1


def test_control_patterns_cover_xenium_terms():
    for term in ("negcontrolprobe", "blank", "antisense", "unassignedcodeword"):
        assert term in mc.CONTROL_PATTERNS


def test_data_driven_qc_keeps_majority_and_drops_controls(tiny_adata, tmp_path):
    n_in = tiny_adata.n_obs
    filtered, summary, qc_obs = mc.data_driven_qc(tiny_adata, outdir=tmp_path)
    assert summary["cells_in"] == n_in
    assert summary["cells_out"] <= n_in
    assert summary["cells_out"] >= 0.5 * n_in            # not over-aggressive
    # The control feature must be dropped from the expression matrix.
    assert not filtered.var_names.str.lower().str.contains("negcontrol").any()
    assert set(qc_obs["status"].unique()) <= {"kept", "removed"}
