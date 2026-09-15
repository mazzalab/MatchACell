"""Test scAdam's boundary with MatchACell without downloading an atlas."""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")

import scparadise_tool as tool  # noqa: E402


@pytest.fixture
def clustered():
    counts = np.array([[1, 4, 2], [3, 1, 5], [2, 3, 1], [5, 1, 2]], dtype="float32")
    a = ad.AnnData(sparse.csr_matrix(np.log1p(counts)))
    a.obs_names = ["c1", "c2", "c3", "c4"]
    a.var_names = ["CD3D", "MS4A1", "LYZ"]
    a.layers["counts"] = sparse.csr_matrix(counts)
    a.layers["lognorm"] = a.X.copy()
    a.obs["leiden_0.5"] = pd.Categorical(["0", "0", "1", "1"])
    a.obs["gradient_score"] = 42.0  # pre-existing metadata must survive
    a.obsm["X_umap"] = np.arange(8).reshape(4, 2).astype(float)
    a.raw = a.copy()
    return a


def test_query_uses_normalized_layer_and_preserves_original(clustered):
    clustered.X = clustered.layers["counts"].copy()
    query, overlap = tool.prepare_query(clustered, ["MS4A1", "CD3D"], "lognorm")
    assert overlap == 1
    assert query.var_names.tolist() == ["CD3D", "MS4A1"]
    np.testing.assert_allclose(query.X.toarray(), np.log1p(clustered.X.toarray()[:, :2]))
    query.X[0, 0] = 100
    assert clustered.layers["lognorm"][0, 0] != 100
    assert clustered.n_vars == 3
    assert not query.layers and not query.obsm and query.raw is None


def test_query_can_use_x(clustered):
    query, _ = tool.prepare_query(clustered, clustered.var_names, "")
    np.testing.assert_array_equal(query.X.toarray(), clustered.X.toarray())


def test_query_rejects_incompatible_data(clustered):
    with pytest.raises(ValueError, match="No genes in common"):
        tool.prepare_query(clustered, ["missing"], "lognorm")
    with pytest.raises(ValueError, match="is missing"):
        tool.prepare_query(clustered, clustered.var_names, "not-a-layer")
    clustered.layers["lognorm"][0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        tool.prepare_query(clustered, clustered.var_names, "lognorm")
    clustered.var_names = ["CD3D", "CD3D", "LYZ"]
    with pytest.raises(ValueError, match="unique"):
        tool.prepare_query(clustered, ["CD3D"], "lognorm")


def test_low_model_overlap_is_reported(clustered):
    with pytest.warns(UserWarning, match="80%"):
        _, overlap = tool.prepare_query(clustered, ["CD3D", "missing"], "lognorm")
    assert overlap == 0.5


@pytest.fixture
def predictions(clustered):
    p = ad.AnnData(obs=pd.DataFrame(index=clustered.obs_names))
    p.obs["scparadise_celltype_l1"] = ["T", "B", "T", "Unk"]
    p.obs["scparadise_celltype_l1_probability"] = [0.9, 0.8, 0.7, 0.99]
    p.obs["scparadise_celltype_l2"] = ["CD4", "B", "CD8", "Unk"]
    p.obs["scparadise_celltype_l2_probability"] = [0.9, 0.4, 0.7, 0.99]
    p.obs["scparadise_unknown"] = [False, False, False, True]
    p.obs["gradient_score"] = [0.1, 0.2, 0.3, 0.9]
    return p


def test_levels_confidence_unknowns_and_mapping(clustered, predictions):
    assert tool.transfer_predictions(clustered, predictions, 0, 0.7) == 2
    assert clustered.obs["cell_type_pred"].tolist() == ["CD4", "Unknown", "CD8", "Unknown"]
    assert clustered.obs["scparadise_celltype_l1"].iloc[-1] == "Unknown"
    assert clustered.obs["scparadise_celltype_l2"].iloc[1] == "B"
    assert (clustered.obs["gradient_score"] == 42).all()
    assert clustered.obs["scparadise_gradient_score"].iloc[-1] == 0.9
    table = tool.cluster_mapping(clustered, "leiden_0.5")
    assert table["cluster"].tolist() == ["0", "1"]
    assert table["Unknown"].tolist() == [1, 1]
    assert table["dominant_pct"].tolist() == [50, 50]
    assert tool.transfer_predictions(clustered, predictions, 1, 0) == 1
    assert clustered.obs["cell_type_pred"].tolist() == ["T", "B", "T", "Unknown"]


def test_prediction_contract_errors(clustered, predictions):
    with pytest.raises(ValueError, match="unavailable"):
        tool.transfer_predictions(clustered, predictions, 3, 0)
    with pytest.raises(ValueError, match="identities or order"):
        tool.transfer_predictions(clustered, predictions[::-1], 0, 0)
    predictions.obs["scparadise_celltype_l2_probability"] = np.nan
    with pytest.raises(ValueError, match="invalid prediction probabilities"):
        tool.transfer_predictions(clustered, predictions, 0, 0)


@pytest.mark.parametrize("kwargs, message", [
    ({"thr": 1.1}, "thr"),
    ({"batch_size": 0}, "positive"),
    ({"n_threads": 0}, "positive"),
    ({"model_dir": ""}, "Configure"),
])
def test_invalid_options_fail_before_loading(kwargs, message):
    options = dict(raw_data_file="missing", verdict_file="missing", model_dir="test_scAdam", output_dir="unused")
    options.update(kwargs)
    with pytest.raises(ValueError, match=message):
        tool.main(**options)


@pytest.mark.parametrize("with_umap", [False, True])
def test_outputs_with_real_scadam_checkpoint(clustered, tmp_path, with_umap):
    """Exercise real serialization/inference; random weights do not test accuracy."""
    scparadise = pytest.importorskip("scparadise")
    torch = pytest.importorskip("torch")
    from sklearn.preprocessing import LabelEncoder

    torch.set_num_threads(1)
    torch.manual_seed(0)
    model = scparadise.scadam.scAdamTransformer(
        gn=3, ncl=[2, 3], ed=8, nc=1, nb=1, nh=2,
        ff_hd=16, classifier_hd=8, dropout=0,
    )
    model.var_names = clustered.var_names.tolist()
    model.celltype_keys = ["coarse", "fine"]
    model.celltype_encoders = {
        "coarse": {"celltype_encoder": LabelEncoder().fit(["Lymphoid", "Myeloid"])},
        "fine": {"celltype_encoder": LabelEncoder().fit(["B", "Mono", "T"])},
    }
    scparadise.scadam.save_model(model, str(tmp_path), "tiny_scAdam", verbose=False)
    if not with_umap:
        del clustered.obsm["X_umap"]
    input_path = tmp_path / "input.h5ad"
    clustered.write_h5ad(input_path)
    verdict = tmp_path / "MatchA_Verdict.txt"
    verdict.write_text("Recommended resolution : leiden_0.5\n")
    out = tmp_path / "annotation output"
    if with_umap:
        completed = subprocess.run(
            [sys.executable, tool.__file__, "--input", str(input_path),
             "--verdict", str(verdict), "--model_dir", str(tmp_path / "tiny_scAdam"),
             "--output", str(out), "--device", "cpu", "--n_threads", "1", "--batch_size", "2"],
            capture_output=True, text=True, timeout=60,
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
    else:
        tool.main(input_path, verdict, str(tmp_path / "tiny_scAdam"), out,
                  device="cpu", n_threads=1, batch_size=2)
    result = ad.read_h5ad(out / "scparadise_annotated.h5ad")
    assert result.obs_names.equals(clustered.obs_names)
    assert result.var_names.equals(clustered.var_names)
    np.testing.assert_array_equal(result.X.toarray(), clustered.X.toarray())
    np.testing.assert_array_equal(result.layers["counts"].toarray(), clustered.layers["counts"].toarray())
    np.testing.assert_array_equal(result.raw.X.toarray(), clustered.raw.X.toarray())
    assert result.uns["scparadise"]["label_level"] == 2
    assert result.obs["prediction_probability"].between(0, 1).all()
    assert set(result.obs["cell_type_pred"].astype(str)) <= {"B", "Mono", "T"}
    table = pd.read_excel(out / "cluster_to_celltype_mapping.xlsx")
    assert len(table) == 2
    count_columns = table.columns.difference(["cluster", "dominant_cell_type", "dominant_pct"])
    assert table[count_columns].to_numpy().sum() == clustered.n_obs
    assert bool(list(out.glob("*.png"))) == with_umap
