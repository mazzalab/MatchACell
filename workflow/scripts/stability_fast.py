"""
MatchACell - fast CPU bootstrap for the Leiden cluster-stability step.

Computes the same estimator as `bootstrap_stability` in
matchacell_cluster_stability.py: each iteration draws a subsample of the cells
(without replacement), rebuilds the kNN graph on the subsample's embedding,
re-runs Leiden, and scores each reference cluster by its best Jaccard overlap
with a subsample cluster. For a given subsample the graph and the partitions
are identical to the per-resolution path (same UMAP connectivities, same
igraph Leiden call and seed); only the work is organised differently:

1. One subsample per iteration is shared by every resolution: the graph is
   built once and Leiden runs at each resolution on it.
2. The subsample kNN is read off one full-data kNN with more neighbours
   (`k_big`): a cell's k nearest neighbours inside the subsample are its first
   k full-data neighbours that survived the draw. With exact kNN (< 4096 cells,
   as in scanpy) this is exactly the subsample's kNN; rows with too few
   survivors are queried directly.
3. Iterations run in parallel worker processes, in batches; after each batch
   the running means are written to `stability/csv/`, so a killed job keeps
   what it computed.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("matchacell")


# =========================================================================== #
# Graph building
# =========================================================================== #
def full_knn(X, k_big, n_jobs=1, seed=0):
    """(indices, distances) of the k_big nearest neighbours of every row, with
    the row itself in column 0 at distance exactly 0."""
    n = X.shape[0]
    k_big = min(k_big, n)
    if n < 4096:  # exact, like scanpy's small-data path
        from sklearn.neighbors import NearestNeighbors

        dist, idx = NearestNeighbors(n_neighbors=k_big, n_jobs=n_jobs).fit(X).kneighbors(X)
    else:
        from pynndescent import NNDescent

        index = NNDescent(X, n_neighbors=k_big, metric="euclidean", random_state=seed,
                          n_jobs=n_jobs, low_memory=True)
        idx, dist = index.neighbor_graph
    idx, dist = idx.astype(np.int64), dist.astype(np.float32)
    # UMAP takes the first non-zero distance as rho, so a round-off self
    # distance (~1e-6) would shift every weight; scanpy zeroes it for the same
    # reason. Rows whose self is missing (duplicate points) get it prepended.
    is_self = idx == np.arange(n)[:, None]
    has_self = is_self.any(axis=1)
    order = np.argsort(~is_self, axis=1, kind="stable")
    idx = np.where(has_self[:, None], np.take_along_axis(idx, order, 1),
                   np.hstack([np.arange(n)[:, None], idx[:, :-1]]))
    dist = np.where(has_self[:, None], np.take_along_axis(dist, order, 1),
                    np.hstack([np.zeros((n, 1), np.float32), dist[:, :-1]]))
    dist[:, 0] = 0.0
    return idx.astype(np.int32), dist


def subsample_knn(full_idx, full_dist, sub, X, k):
    """kNN restricted to `sub` (sorted row ids), in subsample coordinates.
    Returns (indices, distances, number of rows queried directly)."""
    pos = np.full(full_idx.shape[0], -1, dtype=np.int64)
    pos[sub] = np.arange(sub.size)
    cand = pos[full_idx[sub]]                  # -1 where the neighbour was dropped
    keep = cand >= 0
    rank = np.cumsum(keep, axis=1)
    ok = rank[:, -1] >= k
    sel = keep & (rank <= k)
    knn_i = np.empty((sub.size, k), dtype=np.int64)
    knn_d = np.empty((sub.size, k), dtype=np.float32)
    rows = np.flatnonzero(ok)
    knn_i[rows] = cand[rows][sel[rows]].reshape(-1, k)
    knn_d[rows] = full_dist[sub[rows]][sel[rows]].reshape(-1, k)
    bad = np.flatnonzero(~ok)
    if bad.size:
        from sklearn.neighbors import NearestNeighbors

        Xs = X[sub]
        d, i = NearestNeighbors(n_neighbors=k).fit(Xs).kneighbors(Xs[bad])
        knn_i[bad], knn_d[bad] = i, d
        knn_d[bad, 0] = 0.0
    return knn_i, knn_d, int(bad.size)


def connectivities(knn_i, knn_d, k):
    """UMAP fuzzy-simplicial-set weights, as scanpy.pp.neighbors(method='umap')."""
    from scanpy.neighbors import _connectivity

    return _connectivity.umap(knn_i, knn_d, n_obs=knn_i.shape[0], n_neighbors=k)


def leiden_partitions(conn, resolutions, n_iterations=2, seed=0):
    """Leiden at every resolution on one graph; the same call and seeding as
    scanpy.tl.leiden(flavor="igraph", directed=False). Returns codes per
    resolution and the seconds each resolution took."""
    try:
        from scanpy._utils import get_igraph_from_adjacency
        from scanpy._utils.random import set_igraph_random_state
    except ImportError:  # private scanpy helpers moved: go through the public API
        return _leiden_partitions_public(conn, resolutions, n_iterations, seed)

    g = get_igraph_from_adjacency(conn, directed=False)
    out, secs = [], []
    for res in resolutions:
        t0 = time.perf_counter()
        with set_igraph_random_state(seed):
            part = g.community_leiden(objective_function="modularity", weights="weight",
                                      resolution=res, n_iterations=n_iterations)
        out.append(np.asarray(part.membership, dtype=np.int64))
        secs.append(time.perf_counter() - t0)
    return out, secs


def _leiden_partitions_public(conn, resolutions, n_iterations, seed):
    import anndata as ad
    import scanpy as sc

    sub = ad.AnnData(np.zeros((conn.shape[0], 1), dtype="float32"))
    out, secs = [], []
    for res in resolutions:
        t0 = time.perf_counter()
        sc.tl.leiden(sub, resolution=res, key_added="tmp", adjacency=conn, flavor="igraph",
                     n_iterations=n_iterations, directed=False, random_state=seed)
        out.append(sub.obs["tmp"].cat.codes.to_numpy().astype(np.int64))
        secs.append(time.perf_counter() - t0)
    return out, secs


def best_jaccard(ref_codes, n_ref, new_codes):
    """Best Jaccard overlap of each reference cluster with a subsample cluster;
    NaN for reference clusters absent from the subsample."""
    n_new = int(new_codes.max()) + 1
    ct = np.bincount(ref_codes * n_new + new_codes, minlength=n_ref * n_new)
    ct = ct.reshape(n_ref, n_new).astype(float)
    row, col = ct.sum(axis=1), ct.sum(axis=0)
    union = row[:, None] + col[None, :] - ct
    jac = np.divide(ct, union, out=np.zeros_like(ct), where=union > 0)
    best = jac.max(axis=1)
    best[row == 0] = np.nan
    return best


# =========================================================================== #
# Worker processes
# =========================================================================== #
_W = {}


def _init_worker(state):
    try:  # one process per core already: no numba threads on top
        import numba

        numba.set_num_threads(1)
    except Exception:  # pragma: no cover
        pass
    _W.update(state)


def _run_iterations(seeds):
    W = _W
    n, m, k = W["X"].shape[0], W["m"], W["k"]
    out = []
    for seed in seeds:
        sub = np.sort(np.random.default_rng(seed).choice(n, size=m, replace=False))
        knn_i, knn_d, n_direct = subsample_knn(W["full_idx"], W["full_dist"], sub, W["X"], k)
        parts, secs = leiden_partitions(connectivities(knn_i, knn_d, k), W["resolutions"],
                                        W["n_iterations"])
        scores = [best_jaccard(codes[sub], n_ref, lab)
                  for codes, n_ref, lab in zip(W["ref_codes"], W["n_refs"], parts)]
        out.append((scores, secs, n_direct))
    return out


# =========================================================================== #
# Driver
# =========================================================================== #
def _mean_scores(rows):
    """Per-cluster mean over iterations; NaN for a cluster never sampled."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN column
        return np.nanmean(np.vstack(rows), axis=0)


def _write_partial(csv_dir, keys, cats, per_res, n_done):
    for key, c, rows in zip(keys, cats, per_res):
        mean = _mean_scores(rows)
        pd.DataFrame({"Cluster": c, "Jaccard_Stability": mean, "n_iter_done": n_done}) \
            .to_csv(csv_dir / f"stability_{key.replace('leiden_', '')}.csv", index=False)


def bootstrap_stability_all(X, ref_labels, leiden_keys, n_iter, fraction, n_neighbors,
                            seed=0, n_jobs=1, csv_dir: Path | None = None, k_big=None,
                            n_iterations=2):
    """Bootstrap stability for every resolution at once.

    X           : (n, d) embedding the clustering was built on.
    ref_labels  : one length-n label array per entry of `leiden_keys`.
    Returns (scores, leiden_seconds, m): scores[r] maps cluster -> mean Jaccard,
    leiden_seconds[r] is the Leiden time spent on resolution r (summed over
    workers) and m the subsample size.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    n = X.shape[0]
    m = max(int(round(n * fraction)), n_neighbors + 1)
    k_big = k_big or max(3 * n_neighbors, n_neighbors + 15)
    resolutions = [float(k.split("_", 1)[1]) for k in leiden_keys]

    t0 = time.perf_counter()
    full_idx, full_dist = full_knn(X, k_big, n_jobs=n_jobs, seed=seed)
    log.info("Bootstrap: %d-NN of %d cells in %.1fs; %d iterations x %d resolutions "
             "on %d worker(s).", k_big, n, time.perf_counter() - t0, n_iter,
             len(resolutions), n_jobs)

    cats, ref_codes, n_refs = [], [], []
    for lab in ref_labels:
        c, codes = np.unique(np.asarray(lab).astype(str), return_inverse=True)
        cats.append(c)
        ref_codes.append(codes.astype(np.int64))
        n_refs.append(c.size)

    seeds = np.random.default_rng(seed).integers(0, 2**63 - 1, size=n_iter)
    state = dict(X=X, full_idx=full_idx, full_dist=full_dist, ref_codes=ref_codes,
                 n_refs=n_refs, resolutions=resolutions, k=n_neighbors, m=m,
                 n_iterations=n_iterations)
    if csv_dir is not None:
        csv_dir.mkdir(parents=True, exist_ok=True)

    per_res = [[] for _ in resolutions]
    leiden_secs = np.zeros(len(resolutions))
    n_direct = 0
    batch = max(4 * n_jobs, 8)
    done = 0
    # spawn: workers start clean, without the parent's numba/OpenMP threads.
    with ProcessPoolExecutor(max_workers=n_jobs, mp_context=mp.get_context("spawn"),
                             initializer=_init_worker, initargs=(state,)) as ex:
        while done < n_iter:
            todo = seeds[done:done + batch]
            chunks = [todo[i::n_jobs] for i in range(n_jobs) if len(todo[i::n_jobs])]
            for results in ex.map(_run_iterations, chunks):
                for scores, secs, nd in results:
                    for r, s in enumerate(scores):
                        per_res[r].append(s)
                    leiden_secs += secs
                    n_direct += nd
            done += len(todo)
            if csv_dir is not None:
                _write_partial(csv_dir, leiden_keys, cats, per_res, done)
            log.info("  bootstrap: iteration %d/%d (%.0fs)", done, n_iter,
                     time.perf_counter() - t0)
    if n_direct:
        log.info("  %d subsample rows needed a direct kNN query (k_big=%d).",
                 n_direct, k_big)

    scores = []
    for c, rows in zip(cats, per_res):
        mean = _mean_scores(rows)
        scores.append({str(cl): float(v) for cl, v in zip(c, mean)})
    return scores, leiden_secs, m
