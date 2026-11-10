"""Fast CPU bootstrap: same graph and partitions as scanpy, per subsample."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

sc = pytest.importorskip("scanpy")
ad = pytest.importorskip("anndata")
pytest.importorskip("igraph")
import stability_fast as sf  # noqa: E402


@pytest.fixture
def blobs(rng):
    """Three Gaussian blobs in a 10-d embedding (exact kNN regime, < 4096 cells)."""
    centers = rng.normal(0, 6, size=(3, 10))
    X = np.vstack([c + rng.normal(size=(150, 10)) for c in centers]).astype("float32")
    return X


def test_subsample_graph_and_partitions_match_scanpy(blobs, rng):
    X = blobs
    full_idx, full_dist = sf.full_knn(X, 45)
    for _ in range(3):
        sub = np.sort(rng.choice(X.shape[0], size=int(0.8 * X.shape[0]), replace=False))
        ref = ad.AnnData(np.zeros((sub.size, 1), dtype="float32"))
        ref.obsm["X_pca"] = X[sub]
        sc.pp.neighbors(ref, use_rep="X_pca", n_neighbors=15)

        knn_i, knn_d, _ = sf.subsample_knn(full_idx, full_dist, sub, X, 15)
        conn = sf.connectivities(knn_i, knn_d, 15)
        # equal up to float32 round-off of the two kNN distance computations
        assert abs(conn - ref.obsp["connectivities"]).max() < 1e-5

        parts, _ = sf.leiden_partitions(conn, [0.3, 1.0])
        for res, labels in zip([0.3, 1.0], parts):
            sc.tl.leiden(ref, resolution=res, key_added="t", flavor="igraph",
                         n_iterations=2, directed=False)
            # same partition (labels may be numbered differently)
            ct = pd.crosstab(labels, ref.obs["t"].values)
            assert ((ct > 0).sum(axis=1) == 1).all() and ((ct > 0).sum(axis=0) == 1).all()


def test_best_jaccard():
    ref = np.array([0, 0, 0, 1, 1, 1])
    new = np.array([0, 0, 1, 1, 1, 1])
    best = sf.best_jaccard(ref, 3, new)            # cluster 2 is absent
    assert best[0] == pytest.approx(2 / 3)
    assert best[1] == pytest.approx(3 / 4)
    assert np.isnan(best[2])


def test_bootstrap_all_writes_running_csvs(blobs, tmp_path):
    X = blobs
    labels = np.repeat(["0", "1", "2"], 150)
    scores, leiden_secs, m = sf.bootstrap_stability_all(
        X, [labels, labels], ["leiden_0.2", "leiden_0.5"], n_iter=6, fraction=0.8,
        n_neighbors=15, seed=0, n_jobs=2, csv_dir=tmp_path)
    assert m == 360 and len(scores) == 2 and leiden_secs.shape == (2,)
    assert set(scores[0]) == {"0", "1", "2"}
    assert all(s > 0.95 for s in scores[0].values())  # well-separated blobs
    csv = pd.read_csv(tmp_path / "stability_0.2.csv")
    assert list(csv.columns) == ["Cluster", "Jaccard_Stability", "n_iter_done"]
    assert (csv["n_iter_done"] == 6).all()
