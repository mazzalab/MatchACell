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


def _doublet_csv(adata, tmp_path, doublets):
    import pandas as pd
    calls = pd.DataFrame({"cell": adata.obs_names,
                          "doublet_score": [0.9 if c in doublets else 0.05 for c in adata.obs_names],
                          "doublet_class": ["doublet" if c in doublets else "singlet"
                                            for c in adata.obs_names]})
    path = tmp_path / "scDblFinder.csv"
    calls.to_csv(path, index=False)
    return path


def test_doublets_are_a_qc_criterion(tiny_adata, tmp_path):
    path = _doublet_csv(tiny_adata, tmp_path, {"cell3", "cell200"})
    assert mc.read_doublets(tiny_adata, path) == 2
    filtered, summary, qc_obs = mc.data_driven_qc(tiny_adata, outdir=tmp_path)
    assert summary["removed_by_doublet"] == 2
    assert "cell3" not in filtered.obs_names and "cell200" not in filtered.obs_names
    assert (qc_obs.loc[["cell3", "cell200"], "status"] == "removed").all()
    assert {"doublet_score", "doublet_class"} <= set(filtered.obs.columns)


def test_flagged_doublets_stay(tiny_adata, tmp_path):
    path = _doublet_csv(tiny_adata, tmp_path, {"cell3"})
    mc.read_doublets(tiny_adata, path, action="flag")
    filtered, summary, qc_obs = mc.data_driven_qc(tiny_adata, outdir=tmp_path)
    assert summary["removed_by_doublet"] == 0
    assert qc_obs.loc["cell3", "doublet_class"] == "doublet"
    assert not qc_obs["doublet"].any()


def test_doublet_calls_must_cover_every_cell(tiny_adata, tmp_path):
    path = _doublet_csv(tiny_adata[:10].copy(), tmp_path, set())
    with pytest.raises(ValueError, match="no doublet call"):
        mc.read_doublets(tiny_adata, path)
