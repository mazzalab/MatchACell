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

The workers compute UMAP's connectivities and scanpy's igraph Leiden call with
numpy and igraph alone (no umap, numba or scanpy import, which cost ~30 s per
fresh process). tests/test_stability_fast.py checks the graph and partitions
against scanpy; on pbmc3k and PBMC 10k subsamples the weights were bitwise
identical and all 320 partitions the same.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import sys
import time
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np

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


# The workers rebuild the graph with numpy and igraph only. Importing umap and
# JIT-compiling its numba kernels costs ~30 s in every fresh worker process,
# more than the whole bootstrap of a few thousand cells.
def _smooth_knn_dist(dist, k, n_iter=64, tol=1e-5, min_scale=1e-3):
    """umap.umap_.smooth_knn_dist with local_connectivity=1, all rows at once:
    rho is the first non-zero distance, sigma solves
    sum_j exp(-(d_j - rho) / sigma) = log2(k) by bisection."""
    n = dist.shape[0]
    nonzero = dist > 0
    rho = np.where(nonzero.any(axis=1), dist[np.arange(n), nonzero.argmax(axis=1)], 0)
    rho = rho.astype(np.float32)
    shifted = dist[:, 1:] - rho[:, None]                      # float32, as in umap
    above = shifted > 0
    target = np.log2(k)
    lo, hi, mid = np.zeros(n), np.full(n, np.inf), np.ones(n)
    active = np.ones(n, dtype=bool)
    for _ in range(n_iter):
        with np.errstate(over="ignore"):
            psum = np.where(above, np.exp(-(shifted / mid[:, None])), 1.0).sum(axis=1)
        active &= np.abs(psum - target) >= tol
        if not active.any():
            break
        down = active & (psum > target)
        up = active & ~down
        hi[down] = mid[down]
        lo[up] = mid[up]
        mid = np.where(down | (up & np.isfinite(hi)), (lo + hi) / 2,
                       np.where(up, mid * 2, mid))
    sigma = mid.astype(np.float32)
    floor = min_scale * np.where(rho > 0, dist.mean(axis=1), dist.mean())
    return np.maximum(sigma, floor).astype(np.float32), rho


def connectivities(knn_i, knn_d, k):
    """UMAP fuzzy-simplicial-set weights, as scanpy.pp.neighbors(method='umap')
    computes them (umap.fuzzy_simplicial_set, set_op_mix_ratio=1)."""
    from scipy import sparse

    n = knn_i.shape[0]
    knn_d = knn_d.astype(np.float32)
    sigma, rho = _smooth_knn_dist(knn_d, float(k))
    shifted = knn_d - rho[:, None]
    # entries np.where discards may overflow
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        vals = np.where((shifted <= 0) | (sigma[:, None] == 0), np.float32(1),
                        np.exp(-(shifted / sigma[:, None])))
    vals[knn_i == np.arange(n)[:, None]] = 0          # no self-loops
    graph = sparse.coo_matrix((vals.astype(np.float32).ravel(),
                               (np.repeat(np.arange(n), knn_i.shape[1]), knn_i.ravel())),
                              shape=(n, n))
    graph.eliminate_zeros()
    t = graph.transpose()
    graph = graph + t - graph.multiply(t)
    graph.eliminate_zeros()
    return graph.tocsr()


class _IgraphRNG:
    """igraph random number generator over numpy's RandomState, as scanpy
    seeds Leiden (scanpy._utils.random._RNGIgraph)."""

    def __init__(self, seed):
        self._rng = np.random.RandomState(seed)

    def getrandbits(self, k):
        return self._rng.tomaxint() & ((1 << k) - 1)

    def randint(self, a, b):
        return self._rng.randint(a, b + 1)

    def __getattr__(self, attr):
        return getattr(self._rng, "normal" if attr == "gauss" else attr)


def leiden_partitions(conn, resolutions, n_iterations=2, seed=0):
    """Leiden at every resolution on one graph; the same graph, call and
    seeding as scanpy.tl.leiden(flavor="igraph", directed=False). Returns codes
    per resolution and the seconds each resolution took."""
    import random

    import igraph as ig

    # as scanpy.get_igraph_from_adjacency: one edge per non-zero entry
    sources, targets = conn.nonzero()
    g = ig.Graph(directed=False)
    g.add_vertices(conn.shape[0])
    g.add_edges(list(zip(sources, targets)))
    g.es["weight"] = np.asarray(conn[sources, targets]).ravel()
    out, secs = [], []
    for res in resolutions:
        t0 = time.perf_counter()
        ig.set_random_number_generator(_IgraphRNG(seed))
        try:
            part = g.community_leiden(objective_function="modularity", weights="weight",
                                      resolution=res, n_iterations=n_iterations)
        finally:
            ig.set_random_number_generator(random)
        out.append(np.asarray(part.membership, dtype=np.int64))
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
    _W.update(state)


@contextmanager
def _workers(n_jobs, state):
    """A map() over `n_jobs` worker processes, or in this process for one job."""
    if n_jobs == 1:
        _init_worker(state)
        try:
            yield map
        finally:
            _W.clear()
        return
    # spawn: workers start clean, without the parent's numba/OpenMP threads.
    # A spawned worker re-runs the caller's script to rebuild __main__ (for
    # Step 1: scanpy and the plotting libraries, seconds per worker); it only
    # needs this module, so the script's path is hidden while the pool starts
    # its workers (all of them, in Pool()).
    main = sys.modules.get("__main__")
    path = main.__dict__.pop("__file__", None) if main is not None else None
    try:
        pool = mp.get_context("spawn").Pool(n_jobs, _init_worker, (state,))
    finally:
        if path is not None:
            main.__file__ = path
    with pool:
        yield pool.map


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
    import pandas as pd  # parent only: workers skip the import

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
    with _workers(n_jobs, state) as run:
        while done < n_iter:
            todo = seeds[done:done + batch]
            chunks = [todo[i::n_jobs] for i in range(n_jobs) if len(todo[i::n_jobs])]
            for results in run(_run_iterations, chunks):
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
