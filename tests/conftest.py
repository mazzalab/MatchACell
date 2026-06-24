"""Shared pytest fixtures.

The workflow script imports scanpy/anndata at module import time, so test
modules that need it guard with ``pytest.importorskip("scanpy")`` before
importing ``matchacell_cluster_stability``. ``pythonpath`` in pyproject.toml
puts ``workflow/scripts`` on sys.path.
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture
def rng():
    return np.random.default_rng(0)


@pytest.fixture
def tiny_adata(rng):
    """A small two-blob AnnData with raw-ish counts, a couple of MT genes and
    one Xenium-style negative-control feature."""
    ad = pytest.importorskip("anndata")
    n_a, n_b, n_genes = 160, 140, 60
    block_a = rng.poisson(3.0, size=(n_a, n_genes))
    block_b = rng.poisson(6.0, size=(n_b, n_genes))
    X = np.vstack([block_a, block_b]).astype("float32")

    var_names = [f"GENE{i}" for i in range(n_genes)]
    var_names[0] = "MT-CO1"
    var_names[1] = "MT-ND1"
    var_names[2] = "NegControlProbe_00042"
    obs_names = [f"cell{i}" for i in range(n_a + n_b)]

    a = ad.AnnData(X=X)
    a.var_names = var_names
    a.obs_names = obs_names
    return a


@pytest.fixture
def tiny_h5ad(tiny_adata, tmp_path):
    p = tmp_path / "tiny.h5ad"
    tiny_adata.write(p)
    return p
