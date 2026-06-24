#!/usr/bin/env python
"""
MatchACell - Step 1: Data-driven QC + Leiden cluster-stability assessment
=========================================================================

GPU-accelerated (rapids-singlecell) clustering-stability tool with a matcha
theme, rich diagnostics, embedding plots, an interactive Plotly HTML report,
and an automatic "MatchA Verdict" recommending the resolution to annotate.

Pipeline
--------
1.  Load a raw(-ish) .h5ad (post cell-segmentation counts).
2.  Data-driven QC: MAD-based outlier detection (no hard-coded thresholds),
    adaptive to the assay (mito genes and/or Xenium negative controls).
3.  Preprocess on GPU: normalize -> log1p -> (HVG) -> PCA -> neighbours.
4.  Multi-resolution Leiden clustering -> `leiden_<res>` columns.
5.  Embeddings (PCA / t-SNE / UMAP) coloured by every resolution.
6.  Clustree + many QC diagnostic plots.
7.  Bootstrap Jaccard cluster stability.
8.  MatchA Verdict: which resolution to take into annotation, and why.
9.  Interactive Plotly HTML report bundling the key results.

NOTE on theming: the matcha-green theme applies to the CLI and to the HTML
report *chrome* only. The scientific plots keep their default colour maps.

Author: generated for the MatchACell project.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import logging
import os
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib

matplotlib.use("Agg")  # headless / SLURM-safe
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import median_abs_deviation

# Optional dependencies -------------------------------------------------------
try:
    from pyclustree import clustree

    _HAS_CLUSTREE = True
except Exception:  # pragma: no cover
    _HAS_CLUSTREE = False

try:
    import plotly.graph_objects as go
    import plotly.offline as pyo

    _HAS_PLOTLY = True
except Exception:  # pragma: no cover
    _HAS_PLOTLY = False


# =========================================================================== #
# Matcha theme (CLI colours + branding "scent")
# =========================================================================== #
# Truecolor matcha palette (R, G, B).
_MATCHA = (122, 182, 97)      # bright matcha leaf
_DEEP = (74, 124, 60)         # steeped deep green
_LIGHT = (183, 213, 160)      # foam green
_STONE = (140, 150, 120)      # muted sage (for dim text)

MATCHA_SCENT = "Brewed with MatchACell - matcha-grade single-cell consensus"

# Matcha palette used *only* for the HTML report chrome (not the plots).
HTML_THEME = {
    "bg": "#f4f7ee",
    "panel": "#ffffff",
    "ink": "#2f3b24",
    "matcha": "#7AB661",
    "deep": "#4A7C3C",
    "light": "#dcebcd",
    "accent": "#5b8c3e",
    "stone": "#788a66",
}


class Theme:
    """Tiny ANSI helper. Disabled automatically when not a TTY / NO_COLOR."""

    enabled = True

    @classmethod
    def configure(cls, no_color: bool):
        cls.enabled = (
            not no_color
            and os.environ.get("NO_COLOR") is None
            and sys.stdout.isatty()
        )

    @staticmethod
    def _fg(rgb):
        return f"\033[38;2;{rgb[0]};{rgb[1]};{rgb[2]}m"

    @classmethod
    def paint(cls, text, rgb=_MATCHA, bold=False):
        if not cls.enabled:
            return text
        b = "\033[1m" if bold else ""
        return f"{b}{cls._fg(rgb)}{text}\033[0m"

    @classmethod
    def matcha(cls, t, bold=False):
        return cls.paint(t, _MATCHA, bold)

    @classmethod
    def deep(cls, t, bold=False):
        return cls.paint(t, _DEEP, bold)

    @classmethod
    def light(cls, t, bold=False):
        return cls.paint(t, _LIGHT, bold)

    @classmethod
    def dim(cls, t):
        return cls.paint(t, _STONE)


class MatchaFormatter(logging.Formatter):
    """Log formatter that steeps the timestamp/level in matcha green."""

    def format(self, record):
        ts = Theme.dim(datetime.now().strftime("%H:%M:%S"))
        lvl = {
            "INFO": Theme.matcha("leaf"),
            "WARNING": Theme.paint("warn", (200, 170, 60), bold=True),
            "ERROR": Theme.paint("burnt", (200, 90, 60), bold=True),
        }.get(record.levelname, record.levelname.lower())
        return f"{ts} {Theme.deep('|')} {lvl:<5} {Theme.deep('|')} {record.getMessage()}"


def setup_logging():
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(MatchaFormatter())
    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.INFO)


log = logging.getLogger("MatchACell")


BANNER = r"""
 __  __         _          _         _      ____        _  _  _
|  \/  |  __ _ | |_   ___ | |__     / \    / ___|  ___ | || || |
| |\/| | / _` || __| / __|| '_ \   / _ \  | |     / _ \| || || |
| |  | || (_| || |_ | (__ | | | | / ___ \ | |___ |  __/| || ||_|
|_|  |_| \__,_| \__| \___||_| |_|/_/   \_\ \____| \___||_||_|(_)
"""

_TEALEAF = r"""        (
         )
      .-"-.-.
     /  ( )  \      ~ matcha steeped ~
     \   ^   /
      `-...-'
"""


def print_banner():
    print(Theme.matcha(BANNER, bold=True))
    print(Theme.deep("       Multi-Annotator Consensus for single-cell types & states"))
    print(Theme.dim("                   Step 1 - QC & cluster stability"))
    print()


# =========================================================================== #
# Performance profiling
# =========================================================================== #
def gpu_mem_info():
    """Return (used_gb, total_gb) for the active GPU, or (None, None)."""
    try:
        import cupy as cp

        free, total = cp.cuda.runtime.memGetInfo()
        return round((total - free) / 1e9, 3), round(total / 1e9, 3)
    except Exception:  # pragma: no cover
        return None, None


def gpu_name():
    try:
        import cupy as cp

        props = cp.cuda.runtime.getDeviceProperties(0)
        name = props["name"]
        return name.decode() if isinstance(name, bytes) else str(name)
    except Exception:  # pragma: no cover
        return None


class PerfTracker:
    """Lightweight wall-clock profiler for pipeline stages and the bootstrap."""

    def __init__(self, backend_label: str):
        self.backend = backend_label
        self.stages = []          # [{stage, seconds}]
        self.stability_rows = []  # [{resolution, n_clusters, n_iter, seconds, iters_per_s, ...}]
        self.meta = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "backend": backend_label,
            "host": platform.node(),
            "python": platform.python_version(),
            "gpu_name": gpu_name(),
        }
        gused, gtot = gpu_mem_info()
        if gtot is not None:
            self.meta["gpu_total_gb"] = gtot
        self._t_start = time.perf_counter()

    @contextlib.contextmanager
    def stage(self, name: str):
        t = time.perf_counter()
        try:
            yield
        finally:
            secs = round(time.perf_counter() - t, 4)
            self.stages.append({"stage": name, "seconds": secs})
            log.info("%s %s in %s", Theme.dim("[perf]"),
                     name, Theme.matcha(f"{secs:.2f}s"))

    def add_stability(self, resolution, n_clusters, n_iter, seconds, n_cells_sub):
        ips = round(n_iter / seconds, 3) if seconds > 0 else float("nan")
        self.stability_rows.append({
            "resolution": resolution,
            "n_clusters": n_clusters,
            "n_iter": n_iter,
            "seconds": round(seconds, 4),
            "iters_per_s": ips,
            "ms_per_iter": round(1000 * seconds / max(n_iter, 1), 2),
            "n_cells_subsample": int(n_cells_sub),
        })

    def finalize(self, n_cells, n_genes):
        self.meta["total_seconds"] = round(time.perf_counter() - self._t_start, 3)
        self.meta["n_cells"] = int(n_cells)
        self.meta["n_genes"] = int(n_genes)
        gused, _ = gpu_mem_info()
        if gused is not None:
            self.meta["gpu_mem_used_gb_end"] = gused

    def stages_df(self):
        return pd.DataFrame(self.stages)

    def stability_df(self):
        return pd.DataFrame(self.stability_rows)


def write_performance(perf: "PerfTracker", outdir: Path):
    """Write performance CSV/JSON + bar-chart PNGs (default colours)."""
    perfdir = outdir / "performance"
    perfdir.mkdir(parents=True, exist_ok=True)

    stages = perf.stages_df()
    if not stages.empty:
        stages.to_csv(perfdir / "performance_metrics.csv", index=False)
    stab = perf.stability_df()
    if not stab.empty:
        stab.to_csv(perfdir / "stability_timing.csv", index=False)
    (perfdir / "performance_summary.json").write_text(
        json.dumps(perf.meta, indent=2), encoding="utf-8")

    # Stage-duration bar chart
    if not stages.empty:
        fig, ax = plt.subplots(figsize=(8, max(3, 0.5 * len(stages))))
        ax.barh(stages["stage"], stages["seconds"])
        ax.invert_yaxis()
        ax.set_xlabel("seconds")
        ax.set_title(f"Stage wall-clock ({perf.backend})")
        for y, s in enumerate(stages["seconds"]):
            ax.text(s, y, f" {s:.1f}s", va="center")
        fig.tight_layout()
        fig.savefig(perfdir / "stage_durations.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    # Per-resolution stability timing
    if not stab.empty:
        fig, ax = plt.subplots(figsize=(9, 5))
        order = stab.sort_values("resolution", key=lambda s: s.astype(float))
        ax.bar(order["resolution"].astype(str), order["seconds"])
        ax.set_xlabel("Leiden resolution")
        ax.set_ylabel("seconds")
        ax.set_title("Stability bootstrap time per resolution")
        ax2 = ax.twinx()
        ax2.plot(order["resolution"].astype(str), order["iters_per_s"],
                 marker="o", color="black")
        ax2.set_ylabel("iterations / second")
        fig.tight_layout()
        fig.savefig(perfdir / "stability_time_per_resolution.png", dpi=200,
                    bbox_inches="tight")
        plt.close(fig)

    log.info("Performance metrics written to %s", perfdir)


# =========================================================================== #
# Backend selection (GPU rapids-singlecell vs CPU scanpy)
# =========================================================================== #
def init_backend(prefer: str = "auto", pool_fraction: float = 0.5,
                 managed: bool = False):
    """Return (rsc_or_None, use_gpu). prefer in {auto, gpu, cpu}."""
    if prefer == "cpu":
        log.info("Backend: CPU (scanpy), forced by --backend cpu.")
        return None, False

    try:
        import rapids_singlecell as rsc  # noqa: F401

        try:
            import rmm
            from rmm.allocators.cupy import rmm_cupy_allocator
            import cupy as cp

            # A non-managed device pool with an explicit initial size avoids
            # per-allocation cudaMalloc/cudaFree syncs (catastrophic for the
            # thousands of tiny allocations the bootstrap makes) and avoids the
            # page-fault/migration overhead of managed (unified) memory. Sizing
            # the pool up front to a fraction of free VRAM keeps it from growing
            # mid-run. Pass --rmm-managed to opt back into managed memory.
            free, total = cp.cuda.runtime.memGetInfo()
            init_bytes = int(free * float(pool_fraction)) & ~255  # 256B-aligned
            rmm.reinitialize(
                pool_allocator=True,
                managed_memory=bool(managed),
                initial_pool_size=init_bytes,
            )
            cp.cuda.set_allocator(rmm_cupy_allocator)
            log.info("RMM %s pool initialised (%.1f GB of %.1f GB free).",
                     "managed" if managed else "device",
                     init_bytes / 1e9, free / 1e9)
        except Exception as e:  # pragma: no cover
            log.warning("Could not initialise RMM pool (%s); continuing without "
                        "it (expect slower allocation).", e)

        log.info("Backend: GPU (rapids-singlecell). %s", Theme.matcha("[accelerated]"))
        return rsc, True
    except Exception as e:
        if prefer == "gpu":
            log.error("--backend gpu requested but rapids-singlecell import failed: %s", e)
            sys.exit(1)
        log.warning("rapids-singlecell unavailable (%s) -> falling back to CPU scanpy.", e)
        return None, False


# =========================================================================== #
# Data-driven QC
# =========================================================================== #
CONTROL_PATTERNS = [
    "negcontrolprobe",
    "negcontrolcodeword",
    "blank",
    "antisense",
    "unassignedcodeword",
    "deprecatedcodeword",
    "intergenic",
    "genomic",
]


def is_outlier(values: np.ndarray, nmads: float) -> np.ndarray:
    """Boolean MAD-outlier mask: |x - median| > nmads * MAD (both tails)."""
    values = np.asarray(values, dtype=float)
    med = np.median(values)
    mad = median_abs_deviation(values)
    if mad == 0:
        return np.zeros_like(values, dtype=bool)
    return (values < med - nmads * mad) | (values > med + nmads * mad)


def mad_bounds(values: np.ndarray, nmads: float):
    values = np.asarray(values, dtype=float)
    med = np.median(values)
    mad = median_abs_deviation(values)
    return med - nmads * mad, med + nmads * mad


def flag_feature_classes(adata) -> dict:
    up = adata.var_names.str.upper()
    low = adata.var_names.str.lower()
    adata.var["mt"] = up.str.startswith("MT-") | up.str.startswith("MT.")
    pattern = "|".join(CONTROL_PATTERNS)
    adata.var["control"] = low.str.contains(pattern, regex=True)
    n_mt = int(adata.var["mt"].sum())
    n_ctrl = int(adata.var["control"].sum())
    log.info("Feature flags: %d mitochondrial, %d control/probe features.", n_mt, n_ctrl)
    return {"mt": n_mt, "control": n_ctrl}


def data_driven_qc(
    adata,
    nmads: float = 5.0,
    nmads_pct: float = 3.0,
    min_counts_floor: int = 10,
    min_genes_floor: int = 5,
    mt_pct_hard: float = 20.0,
    outdir: Path | None = None,
):
    """MAD-based adaptive QC. Returns (filtered_adata, summary, qc_obs)."""
    flags = flag_feature_classes(adata)
    qc_vars = [v for v in ("mt", "control") if flags[v] > 0]
    percent_top = [20] if adata.n_vars >= 20 else None
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=qc_vars, percent_top=percent_top, log1p=True, inplace=True
    )

    n0 = adata.n_obs
    top_col = next((c for c in adata.obs.columns if c.startswith("pct_counts_in_top_")), None)

    adata.obs["outlier_counts"] = is_outlier(adata.obs["log1p_total_counts"], nmads)
    adata.obs["outlier_genes"] = is_outlier(adata.obs["log1p_n_genes_by_counts"], nmads)
    adata.obs["outlier_top"] = (
        is_outlier(adata.obs[top_col], nmads) if top_col else False
    )
    adata.obs["outlier"] = (
        adata.obs["outlier_counts"] | adata.obs["outlier_genes"] | adata.obs["outlier_top"]
    )

    if "pct_counts_mt" in adata.obs:
        adata.obs["mt_outlier"] = is_outlier(adata.obs["pct_counts_mt"], nmads_pct) | (
            adata.obs["pct_counts_mt"] > mt_pct_hard
        )
    else:
        adata.obs["mt_outlier"] = False

    if "pct_counts_control" in adata.obs:
        adata.obs["control_outlier"] = is_outlier(adata.obs["pct_counts_control"], nmads_pct)
    else:
        adata.obs["control_outlier"] = False

    adata.obs["below_floor"] = (
        (adata.obs["total_counts"] < min_counts_floor)
        | (adata.obs["n_genes_by_counts"] < min_genes_floor)
    ).values

    keep = (
        ~adata.obs["outlier"]
        & ~adata.obs["mt_outlier"]
        & ~adata.obs["control_outlier"]
        & ~adata.obs["below_floor"]
    )

    qc_obs = adata.obs.copy()
    qc_obs["status"] = np.where(keep.values, "kept", "removed")

    if outdir is not None:
        qc_diagnostics(adata, keep, nmads, nmads_pct, outdir)

    adata = adata[keep.values].copy()

    if adata.var["control"].any():
        n_ctrl = int(adata.var["control"].sum())
        adata = adata[:, ~adata.var["control"].values].copy()
        log.info("Dropped %d control features from the expression matrix.", n_ctrl)

    n1 = adata.n_obs
    summary = {
        "cells_in": n0,
        "cells_out": n1,
        "removed": n0 - n1,
        "pct_removed": round(100 * (n0 - n1) / max(n0, 1), 2),
        "removed_by_counts": int(qc_obs["outlier_counts"].sum()),
        "removed_by_genes": int(qc_obs["outlier_genes"].sum()),
        "removed_by_mt": int(qc_obs["mt_outlier"].sum()),
        "removed_by_control": int(qc_obs["control_outlier"].sum()),
        "removed_by_floor": int(qc_obs["below_floor"].sum()),
    }
    log.info(
        "QC: kept %s / %d cells (removed %d, %.2f%%).",
        Theme.matcha(str(n1), bold=True),
        n0,
        n0 - n1,
        summary["pct_removed"],
    )
    return adata, summary, qc_obs


# =========================================================================== #
# QC diagnostic plots (many; default colours)
# =========================================================================== #
def qc_diagnostics(adata, keep_mask, nmads, nmads_pct, outdir: Path):
    """A battery of QC diagnostics. Plots use default matplotlib colours."""
    qcdir = outdir / "qc"
    qcdir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs.copy()
    obs["status"] = np.where(keep_mask.values, "kept", "removed")

    metrics = [m for m in ["total_counts", "n_genes_by_counts", "pct_counts_mt",
                           "pct_counts_control"] if m in obs.columns]
    top_col = next((c for c in obs.columns if c.startswith("pct_counts_in_top_")), None)
    if top_col:
        metrics.append(top_col)

    # 1) Violin of each metric, kept vs removed
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
    axes = np.atleast_1d(axes)
    for ax, m in zip(axes, metrics):
        sns.violinplot(data=obs, x="status", y=m, hue="status", legend=False, ax=ax, cut=0)
        ax.set_title(m)
    fig.suptitle("QC metrics: kept vs removed")
    fig.tight_layout()
    fig.savefig(qcdir / "qc_violin_kept_vs_removed.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 2) Histograms with MAD threshold lines
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
    axes = np.atleast_1d(axes)
    for ax, m in zip(axes, metrics):
        logscale = m in ("total_counts", "n_genes_by_counts")
        vals = np.log1p(obs[m]) if logscale else obs[m]
        sns.histplot(vals, bins=80, ax=ax)
        nm = nmads_pct if m.startswith("pct_counts_") else nmads
        lo, hi = mad_bounds(vals, nm)
        for b in (lo, hi):
            ax.axvline(b, color="black", linestyle="--", linewidth=1)
        ax.set_title(("log1p " if logscale else "") + m)
    fig.suptitle(f"Distributions with MAD bounds (nmads={nmads}/{nmads_pct})")
    fig.tight_layout()
    fig.savefig(qcdir / "qc_hist_mad_bounds.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 3) Counts vs genes scatter, coloured by status and by mito%/control%
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.scatterplot(data=obs, x="total_counts", y="n_genes_by_counts", hue="status",
                    s=6, alpha=0.5, ax=axes[0])
    axes[0].set(xscale="log", yscale="log", title="Counts vs genes (QC status)")
    color_metric = "pct_counts_mt" if "pct_counts_mt" in obs else (
        "pct_counts_control" if "pct_counts_control" in obs else None)
    if color_metric:
        sca = axes[1].scatter(obs["total_counts"], obs["n_genes_by_counts"],
                              c=obs[color_metric], s=6, alpha=0.6)
        axes[1].set(xscale="log", yscale="log", title=f"coloured by {color_metric}")
        fig.colorbar(sca, ax=axes[1], label=color_metric)
    else:
        axes[1].axis("off")
    fig.tight_layout()
    fig.savefig(qcdir / "qc_counts_vs_genes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 4) Barcode-rank / knee plot
    fig, ax = plt.subplots(figsize=(6, 5))
    ranked = np.sort(obs["total_counts"].values)[::-1]
    ax.plot(np.arange(1, len(ranked) + 1), ranked)
    ax.set(xscale="log", yscale="log", xlabel="cell rank",
           ylabel="total counts", title="Barcode-rank (knee) plot")
    fig.tight_layout()
    fig.savefig(qcdir / "qc_barcode_rank.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 5) Per-filter removal bar chart
    cats = {
        "low/high counts": "outlier_counts",
        "low/high genes": "outlier_genes",
        "top-genes %": "outlier_top",
        "mito %": "mt_outlier",
        "control %": "control_outlier",
        "below floor": "below_floor",
    }
    counts = {k: int(obs[v].sum()) for k, v in cats.items()
              if v in obs.columns and obs[v].dtype != object}
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(list(counts.keys()), list(counts.values()))
    ax.set_ylabel("cells flagged")
    ax.set_title("Cells flagged per QC criterion (criteria can overlap)")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(qcdir / "qc_removal_by_criterion.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 6) Highest-expressed genes (scanpy)
    try:
        sc.pl.highest_expr_genes(adata, n_top=20, show=False)
        plt.savefig(qcdir / "qc_highest_expr_genes.png", dpi=200, bbox_inches="tight")
        plt.close()
    except Exception as e:  # pragma: no cover
        log.warning("highest_expr_genes plot skipped: %s", e)

    log.info("QC diagnostics written to %s", qcdir)


# =========================================================================== #
# Preprocessing + clustering + embeddings
# =========================================================================== #
def preprocess_and_cluster(adata, rsc, use_gpu, resolutions, n_pcs, n_neighbors,
                           use_hvg, do_tsne, do_umap, target_sum=None):
    sc.pp.filter_genes(adata, min_cells=1)
    adata.layers["counts"] = adata.X.copy()

    if use_gpu:
        rsc.get.anndata_to_GPU(adata)

    pp = rsc.pp if use_gpu else sc.pp
    tl = rsc.tl if use_gpu else sc.tl

    pp.normalize_total(adata, target_sum=target_sum)
    pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    if use_hvg:
        pp.highly_variable_genes(adata, n_top_genes=2000)
        # Subset to HVG for PCA/clustering. Note: this also subsets the `counts`
        # and `lognorm` layers to the HVG set. We deliberately do NOT set
        # `adata.raw` here because raw would hold GPU arrays and break the .h5ad
        # write; full-gene matrices are not needed for the stability step.
        adata = adata[:, adata.var["highly_variable"]].copy()

    pp.pca(adata, n_comps=n_pcs)
    pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    for res in resolutions:
        key = f"leiden_{res}"
        if use_gpu:
            tl.leiden(adata, resolution=float(res), key_added=key)
        else:
            tl.leiden(adata, resolution=float(res), key_added=key,
                      flavor="igraph", n_iterations=2, directed=False)
        log.info("Leiden res=%s -> %d clusters.", res, adata.obs[key].nunique())

    # Embeddings (best-effort; a failure just skips that embedding)
    if do_umap:
        try:
            tl.umap(adata)
            log.info("UMAP computed.")
        except Exception as e:  # pragma: no cover
            log.warning("UMAP skipped: %s", e)
    if do_tsne:
        try:
            tl.tsne(adata)
            log.info("t-SNE computed.")
        except Exception as e:  # pragma: no cover
            log.warning("t-SNE skipped: %s", e)

    if use_gpu:
        rsc.get.anndata_to_CPU(adata)
    return adata


def plot_embeddings(adata, leiden_keys, outdir: Path):
    """PCA / t-SNE / UMAP coloured by every resolution. Default scanpy colours."""
    embdir = outdir / "embeddings"
    embdir.mkdir(parents=True, exist_ok=True)

    bases = [(name, "X_" + name) for name in ("pca", "tsne", "umap")]
    available = [(n, k) for n, k in bases if k in adata.obsm]
    if not available:
        log.warning("No embeddings available to plot.")
        return

    for res_key in leiden_keys:
        n = len(available)
        fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
        axes = np.atleast_1d(axes)
        for ax, (name, _) in zip(axes, available):
            sc.pl.embedding(adata, basis=name, color=res_key, ax=ax, show=False,
                            frameon=False, legend_loc="right margin", title=name.upper())
        fig.suptitle(res_key)
        fig.tight_layout()
        fig.savefig(embdir / f"embeddings_{res_key}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    # QC-coloured embeddings on the richest available basis (prefer UMAP)
    qc_color = [m for m in ["total_counts", "n_genes_by_counts", "pct_counts_mt",
                            "pct_counts_control"] if m in adata.obs]
    if qc_color:
        basis = available[-1][0]
        fig, axes = plt.subplots(1, len(qc_color), figsize=(5 * len(qc_color), 4.5))
        axes = np.atleast_1d(axes)
        for ax, m in zip(axes, qc_color):
            sc.pl.embedding(adata, basis=basis, color=m, ax=ax, show=False,
                            frameon=False, title=m)
        fig.suptitle(f"QC metrics on {basis.upper()}")
        fig.tight_layout()
        fig.savefig(embdir / "embeddings_qc_metrics.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    log.info("Embedding plots written to %s", embdir)


# =========================================================================== #
# Bootstrap Jaccard stability
# =========================================================================== #
def _cluster_subsample(pca_sub, resolution, rsc, use_gpu, n_neighbors):
    """Re-cluster one subsample. `pca_sub` is already a device (cupy) array when
    use_gpu, else a host numpy array - no extra host<->device copy here."""
    sub = ad.AnnData(np.zeros((pca_sub.shape[0], 1), dtype="float32"))
    sub.obsm["X_pca"] = pca_sub
    if use_gpu:
        rsc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=n_neighbors)
        rsc.tl.leiden(sub, resolution=float(resolution), key_added="tmp")
    else:
        sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=n_neighbors)
        sc.tl.leiden(sub, resolution=float(resolution), key_added="tmp",
                     flavor="igraph", n_iterations=2, directed=False)
    return sub.obs["tmp"].astype(str).values


def bootstrap_stability(adata, res_key, resolution, rsc, use_gpu, n_iter,
                        fraction, n_neighbors, rng, X_dev=None):
    """Bootstrap Jaccard stability for one resolution.

    When `X_dev` is provided (a cupy array of the full PCA matrix already on the
    GPU), each iteration gathers its subsample on-device (`X_dev[idx]`) instead
    of copying a fresh PCA slice from host to device every time. This collapses
    ~n_iter host->device transfers into the single up-front transfer that built
    `X_dev`.
    """
    labels = adata.obs[res_key].astype(str).values
    clusters = pd.unique(labels)
    scores = {c: [] for c in clusters}

    n = adata.n_obs
    k = max(int(round(n * fraction)), n_neighbors + 1)

    use_dev = use_gpu and X_dev is not None
    if use_dev:
        import cupy as cp
    else:
        X_host = np.asarray(adata.obsm["X_pca"], dtype="float32")

    for it in range(n_iter):
        idx = rng.choice(n, size=k, replace=False)
        orig_sub = labels[idx]
        if use_dev:
            pca_sub = X_dev[cp.asarray(idx)]          # gather on device
        else:
            pca_sub = X_host[idx]
        temp = _cluster_subsample(pca_sub, resolution, rsc, use_gpu, n_neighbors)

        ct = pd.crosstab(pd.Series(orig_sub, name="orig"),
                         pd.Series(temp, name="temp"))
        inter = ct.values.astype(float)
        row = ct.sum(axis=1).to_numpy()[:, None]
        col = ct.sum(axis=0).to_numpy()[None, :]
        union = row + col - inter
        with np.errstate(divide="ignore", invalid="ignore"):
            jac = np.where(union > 0, inter / union, 0.0)
        best = jac.max(axis=1)

        for c, b in zip(ct.index.tolist(), best):
            scores[c].append(float(b))

        if (it + 1) % max(1, n_iter // 5) == 0:
            log.info("  %s: iteration %d/%d", res_key, it + 1, n_iter)

    return {c: (float(np.mean(v)) if v else np.nan) for c, v in scores.items()}, k


def run_all_stability(adata, leiden_keys, rsc, use_gpu, n_iter, fraction,
                      n_neighbors, outdir, rng, perf=None):
    stabdir = outdir / "stability"
    stabdir.mkdir(parents=True, exist_ok=True)
    excel_path = stabdir / "Cluster_Stability_Summary.xlsx"

    # Move the full PCA matrix to the GPU ONCE and reuse it for every resolution
    # and every bootstrap iteration (the PCA is identical across resolutions).
    X_dev = None
    if use_gpu:
        try:
            import cupy as cp

            X_dev = cp.asarray(np.asarray(adata.obsm["X_pca"], dtype="float32"))
            log.info("PCA matrix moved to GPU once (%s rows) for the bootstrap.",
                     X_dev.shape[0])
            # Warm-up: compile/init cuML+cuGraph kernels so the first timed
            # resolution is not charged for one-off JIT/allocation costs.
            try:
                _ = _cluster_subsample(
                    X_dev[cp.asarray(rng.choice(adata.n_obs,
                          size=min(adata.n_obs, max(50, n_neighbors + 1)),
                          replace=False))],
                    float(leiden_keys[0].split("_", 1)[1]), rsc, True, n_neighbors)
                log.info("GPU kernels warmed up.")
            except Exception as e:  # pragma: no cover
                log.warning("Warm-up skipped: %s", e)
        except Exception as e:  # pragma: no cover
            log.warning("Could not stage PCA on GPU (%s); per-iter transfer used.", e)
            X_dev = None

    rows = []
    with pd.ExcelWriter(excel_path) as writer:
        for res_key in leiden_keys:
            resolution = float(res_key.split("_", 1)[1])
            log.info("Stability for %s ...", res_key)
            t0 = time.perf_counter()
            res, k = bootstrap_stability(adata, res_key, resolution, rsc, use_gpu,
                                         n_iter, fraction, n_neighbors, rng, X_dev=X_dev)
            secs = time.perf_counter() - t0
            if perf is not None:
                perf.add_stability(res_key.replace("leiden_", ""),
                                   len(res), n_iter, secs, k)
            log.info("  %s done in %s (%.1f iters/s).", res_key,
                     Theme.matcha(f"{secs:.1f}s"), n_iter / max(secs, 1e-9))
            df = pd.DataFrame.from_dict(res, orient="index", columns=["Jaccard_Stability"])
            df.index.name = "Cluster"
            df.to_excel(writer, sheet_name=res_key[:31])
            for cluster, score in res.items():
                rows.append({"Resolution": res_key.replace("leiden_", ""),
                             "Cluster": cluster, "Jaccard_Stability": score})

    # Free the device matrix promptly.
    if X_dev is not None:
        try:
            del X_dev
            import cupy as cp

            cp.get_default_memory_pool().free_all_blocks()
        except Exception:  # pragma: no cover
            pass

    master = pd.DataFrame(rows)

    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    order = sorted(master["Resolution"].unique(), key=float)
    sns.boxplot(data=master, x="Resolution", y="Jaccard_Stability", order=order,
                hue="Resolution", legend=False, showfliers=False)
    sns.stripplot(data=master, x="Resolution", y="Jaccard_Stability", order=order,
                  color="black", size=4, jitter=True, alpha=0.6)
    plt.axhline(0.85, color="forestgreen", linestyle="--", linewidth=1.5,
                label="High stability (>0.85)")
    plt.axhline(0.60, color="firebrick", linestyle="--", linewidth=1.5,
                label="Risk threshold (<0.60)")
    plt.title("Clustering stability (bootstrap Jaccard index)")
    plt.ylabel("Jaccard stability index")
    plt.xlabel("Leiden resolution")
    plt.ylim(0, 1.05)
    plt.legend(loc="lower left")
    plt.savefig(stabdir / "global_stability_comparison_boxplot.png", dpi=300,
                bbox_inches="tight")
    plt.close()
    return master


# =========================================================================== #
# MatchA Verdict
# =========================================================================== #
def resolution_summary(master: pd.DataFrame) -> pd.DataFrame:
    g = master.groupby("Resolution")["Jaccard_Stability"]
    summary = pd.DataFrame({
        "n_clusters": g.size(),
        "median": g.median(),
        "mean": g.mean(),
        "min": g.min(),
        "frac_high": g.apply(lambda s: float((s >= 0.85).mean())),
        "frac_risk": g.apply(lambda s: float((s < 0.60).mean())),
    })
    summary = summary.loc[sorted(summary.index, key=float)]
    return summary


def compute_verdict(master: pd.DataFrame, high=0.85, risk=0.60,
                    min_median=0.75, max_risk_frac=0.25) -> dict:
    """
    Recommend a resolution for annotation.

    Rule: among resolutions whose median stability >= `min_median` and whose
    fraction of risky clusters (< `risk`) <= `max_risk_frac`, pick the FINEST
    (most clusters) - finer granularity gives more cell-type/state resolution
    for annotation, provided it stays stable. If none qualify, fall back to the
    resolution with the highest median stability.
    """
    summary = resolution_summary(master)

    acceptable = summary[(summary["median"] >= min_median)
                         & (summary["frac_risk"] <= max_risk_frac)]
    if not acceptable.empty:
        rec = acceptable.sort_values(["n_clusters", "median"],
                                     ascending=[False, False]).index[0]
        basis = "finest stable resolution"
    else:
        rec = summary["median"].idxmax()
        basis = "fallback: highest median stability (no resolution met the bar)"

    rec_row = summary.loc[rec]
    unstable = (
        master[(master["Resolution"] == rec) & (master["Jaccard_Stability"] < risk)]
        .sort_values("Jaccard_Stability")[["Cluster", "Jaccard_Stability"]]
    )

    return {
        "summary": summary,
        "recommended": rec,
        "basis": basis,
        "n_clusters": int(rec_row["n_clusters"]),
        "median": float(rec_row["median"]),
        "frac_high": float(rec_row["frac_high"]),
        "frac_risk": float(rec_row["frac_risk"]),
        "unstable_clusters": unstable,
        "thresholds": {"high": high, "risk": risk,
                       "min_median": min_median, "max_risk_frac": max_risk_frac},
    }


def render_verdict_text(v: dict) -> str:
    lines = []
    lines.append("MatchA Verdict")
    lines.append("=" * 60)
    lines.append(f"Recommended resolution : leiden_{v['recommended']}")
    lines.append(f"Number of clusters     : {v['n_clusters']}")
    lines.append(f"Median Jaccard         : {v['median']:.3f}")
    lines.append(f"Clusters >= {v['thresholds']['high']:.2f}      : {100*v['frac_high']:.0f}%")
    lines.append(f"Clusters <  {v['thresholds']['risk']:.2f}      : {100*v['frac_risk']:.0f}%")
    lines.append(f"Selection basis        : {v['basis']}")
    lines.append("")
    if len(v["unstable_clusters"]):
        lines.append("Watch these clusters during annotation (below risk threshold):")
        for _, r in v["unstable_clusters"].iterrows():
            lines.append(f"  - cluster {r['Cluster']}: Jaccard {r['Jaccard_Stability']:.3f}")
        lines.append("  (consider merging, re-examining markers, or sub-clustering.)")
    else:
        lines.append("No clusters fall below the risk threshold at this resolution.")
    lines.append("")
    lines.append("Per-resolution summary:")
    lines.append(v["summary"].round(3).to_string())
    lines.append("")
    lines.append(MATCHA_SCENT)
    return "\n".join(lines)


def print_verdict(v: dict):
    title = f"  MatchA Verdict: annotate at leiden_{v['recommended']}  "
    bar = Theme.matcha("=" * (len(title) + 2), bold=True)
    print()
    print(bar)
    print(Theme.matcha("| ", bold=True) + Theme.deep(title, bold=True)
          + Theme.matcha(" |", bold=True))
    print(bar)
    print(Theme.deep(f"  clusters: {v['n_clusters']}   "
                     f"median Jaccard: {v['median']:.3f}   "
                     f">=high: {100*v['frac_high']:.0f}%   "
                     f"<risk: {100*v['frac_risk']:.0f}%"))
    print(Theme.dim(f"  basis: {v['basis']}"))
    if len(v["unstable_clusters"]):
        watch = ", ".join(str(c) for c in v["unstable_clusters"]["Cluster"].tolist())
        print(Theme.paint(f"  watch clusters: {watch}", (200, 170, 60)))
    print()


# =========================================================================== #
# Transitions
# =========================================================================== #
def run_transitions(adata, leiden_keys, outdir):
    trdir = outdir / "transitions"
    trdir.mkdir(parents=True, exist_ok=True)
    all_transitions = {}
    maxima = []
    for i in range(len(leiden_keys) - 1):
        ct = pd.crosstab(adata.obs[leiden_keys[i]], adata.obs[leiden_keys[i + 1]])
        maxima.append(ct.values.max())
    vmax = max(maxima) if maxima else 1

    for i in range(len(leiden_keys) - 1):
        parent, child = leiden_keys[i], leiden_keys[i + 1]
        log.info("Transition: %s -> %s", parent, child)
        ct = pd.crosstab(adata.obs[parent], adata.obs[child])
        all_transitions[f"{parent}_to_{child}"] = ct

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(ct.values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)
        ax.set_xticks(np.arange(len(ct.columns)))
        ax.set_yticks(np.arange(len(ct.index)))
        ax.set_xticklabels(ct.columns, rotation=90)
        ax.set_yticklabels(ct.index)
        ax.set_xlabel(child)
        ax.set_ylabel(parent)
        ax.set_title(f"{parent} -> {child}")
        cbar = fig.colorbar(im)
        cbar.set_label("Cell count")
        fig.tight_layout()
        fig.savefig(trdir / f"{parent}_to_{child}_heatmap_counts.png", dpi=300)
        plt.close(fig)

    excel_path = trdir / "All_cluster_transitions.xlsx"
    with pd.ExcelWriter(excel_path) as writer:
        for name, df in all_transitions.items():
            df.to_excel(writer, sheet_name=name[:30])
    return all_transitions


# =========================================================================== #
# Interactive Plotly HTML report
# =========================================================================== #
def _fig_div(fig):
    return fig.to_html(full_html=False, include_plotlyjs=False,
                       default_width="100%", default_height="540px")


def build_html_report(adata, master, verdict, qc_obs, summary, transitions,
                      outdir: Path, leiden_keys, perf=None):
    if not _HAS_PLOTLY:
        log.warning("plotly not installed -> skipping HTML report. `pip install plotly`.")
        return None

    divs = []

    # --- QC scatter (interactive) ---
    if qc_obs is not None and "total_counts" in qc_obs:
        cmetric = "pct_counts_mt" if "pct_counts_mt" in qc_obs else (
            "pct_counts_control" if "pct_counts_control" in qc_obs else "n_genes_by_counts")
        f = go.Figure(go.Scattergl(
            x=qc_obs["total_counts"], y=qc_obs["n_genes_by_counts"], mode="markers",
            marker=dict(size=4, color=qc_obs[cmetric], colorscale="Viridis",
                        showscale=True, colorbar=dict(title=cmetric), opacity=0.6),
            text=qc_obs["status"],
            hovertemplate="counts=%{x}<br>genes=%{y}<br>%{text}<extra></extra>"))
        f.update_layout(title="QC: counts vs genes", xaxis_type="log", yaxis_type="log",
                        xaxis_title="total_counts", yaxis_title="n_genes_by_counts")
        divs.append(("Quality control", _fig_div(f)))

    # --- Stability box (interactive) ---
    if master is not None and len(master):
        order = sorted(master["Resolution"].unique(), key=float)
        f = go.Figure()
        for r in order:
            f.add_trace(go.Box(y=master.loc[master["Resolution"] == r, "Jaccard_Stability"],
                               name=r, boxpoints="all", jitter=0.4, pointpos=0))
        f.add_hline(y=0.85, line_dash="dash", line_color="green",
                    annotation_text="high 0.85")
        f.add_hline(y=0.60, line_dash="dash", line_color="firebrick",
                    annotation_text="risk 0.60")
        f.update_layout(title="Cluster stability by resolution", showlegend=False,
                        yaxis_title="Jaccard", xaxis_title="resolution")
        divs.append(("Stability", _fig_div(f)))

    # --- Interactive embedding with a resolution switcher ---
    basis = next((b for b in ("X_umap", "X_tsne", "X_pca") if b in adata.obsm), None)
    if basis is not None and leiden_keys:
        xy = np.asarray(adata.obsm[basis])
        f = go.Figure()
        rec_key = f"leiden_{verdict['recommended']}" if verdict else leiden_keys[0]
        for key in leiden_keys:
            codes = adata.obs[key].astype("category").cat.codes.values
            f.add_trace(go.Scattergl(
                x=xy[:, 0], y=xy[:, 1], mode="markers", name=key,
                visible=(key == rec_key),
                marker=dict(size=4, color=codes, colorscale="Turbo", showscale=False),
                text=adata.obs[key].astype(str),
                hovertemplate="cluster %{text}<extra></extra>"))
        buttons = []
        for i, key in enumerate(leiden_keys):
            vis = [j == i for j in range(len(leiden_keys))]
            buttons.append(dict(label=key.replace("leiden_", "res "), method="update",
                                args=[{"visible": vis},
                                      {"title": f"{basis[2:].upper()} - {key}"}]))
        f.update_layout(
            title=f"{basis[2:].upper()} - {rec_key}",
            updatemenus=[dict(buttons=buttons, direction="down", x=1.0, xanchor="right",
                              y=1.15, yanchor="top", showactive=True)],
            xaxis_title=f"{basis[2:].upper()}1", yaxis_title=f"{basis[2:].upper()}2")
        divs.append(("Embedding (use the menu, top-right, to switch resolution)", _fig_div(f)))

    # --- Transition heatmap for the recommended boundary ---
    if transitions:
        rec = verdict["recommended"] if verdict else None
        pick = None
        if rec is not None:
            pick = next((n for n in transitions if n.startswith(f"leiden_{rec}_to_")
                         or n.endswith(f"_to_leiden_{rec}")), None)
        name = pick or next(iter(transitions))
        ct = transitions[name]
        f = go.Figure(go.Heatmap(z=ct.values, x=[str(c) for c in ct.columns],
                                 y=[str(r) for r in ct.index], colorscale="YlGnBu"))
        f.update_layout(title=f"Transition: {name}",
                        xaxis_title=name.split("_to_")[-1],
                        yaxis_title=name.split("_to_")[0])
        divs.append(("Cluster transitions", _fig_div(f)))

    # --- Performance (stage durations + per-resolution throughput) ---
    if perf is not None:
        stages = perf.stages_df()
        if not stages.empty:
            f = go.Figure(go.Bar(x=stages["seconds"], y=stages["stage"],
                                 orientation="h"))
            f.update_layout(title=f"Stage wall-clock ({perf.backend})",
                            xaxis_title="seconds",
                            yaxis=dict(autorange="reversed"))
            divs.append(("Performance - stage durations", _fig_div(f)))
        stab = perf.stability_df()
        if not stab.empty:
            order = stab.sort_values("resolution", key=lambda s: s.astype(float))
            f = go.Figure()
            f.add_trace(go.Bar(x=order["resolution"].astype(str),
                               y=order["seconds"], name="seconds"))
            f.add_trace(go.Scatter(x=order["resolution"].astype(str),
                                   y=order["iters_per_s"], name="iters/s",
                                   mode="lines+markers", yaxis="y2"))
            f.update_layout(
                title="Stability bootstrap time per resolution",
                xaxis_title="resolution", yaxis_title="seconds",
                yaxis2=dict(title="iterations / second", overlaying="y", side="right"),
                legend=dict(orientation="h"))
            divs.append(("Performance - stability throughput", _fig_div(f)))

    html = _assemble_html(divs, verdict, summary, perf)
    out = outdir / "MatchACell_report.html"
    out.write_text(html, encoding="utf-8")
    log.info("Interactive report written to %s", out)
    return out


def _assemble_html(divs, verdict, summary, perf=None):
    T = HTML_THEME
    plotly_js = pyo.get_plotlyjs()

    perf_strip = ""
    if perf is not None:
        m = perf.meta
        bits = [f"backend: {m.get('backend', '?')}"]
        if m.get("gpu_name"):
            bits.append(f"GPU: {m['gpu_name']}")
        if m.get("n_cells") is not None:
            bits.append(f"{m['n_cells']:,} cells x {m.get('n_genes', '?')} genes")
        if m.get("total_seconds") is not None:
            bits.append(f"total {m['total_seconds']:.1f}s")
        perf_strip = '<div class="scent">' + " &middot; ".join(bits) + "</div>"

    verdict_html = ""
    if verdict is not None:
        watch = ", ".join(str(c) for c in verdict["unstable_clusters"]["Cluster"].tolist()) \
            if len(verdict["unstable_clusters"]) else "none"
        verdict_html = f"""
        <div class="verdict">
          <div class="verdict-title">MatchA Verdict</div>
          <div class="verdict-main">Annotate at <b>leiden_{verdict['recommended']}</b></div>
          <div class="verdict-sub">
            {verdict['n_clusters']} clusters &middot; median Jaccard {verdict['median']:.3f}
            &middot; {100*verdict['frac_high']:.0f}% high &middot; {100*verdict['frac_risk']:.0f}% risky
          </div>
          <div class="verdict-basis">{verdict['basis']}</div>
          <div class="verdict-watch">watch clusters: {watch}</div>
        </div>"""

    summary_html = summary.round(3).to_html(classes="summary", border=0) \
        if summary is not None else ""

    sections = "\n".join(
        f'<section><h2>{title}</h2><div class="card">{body}</div></section>'
        for title, body in divs
    )

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MatchACell report</title>
<script>{plotly_js}</script>
<style>
  body {{ margin:0; background:{T['bg']}; color:{T['ink']};
         font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif; }}
  header {{ background:linear-gradient(135deg,{T['matcha']},{T['deep']});
           color:#fff; padding:28px 32px; }}
  header h1 {{ margin:0; font-size:26px; letter-spacing:.5px; }}
  header .scent {{ opacity:.92; margin-top:6px; font-size:13px; }}
  main {{ max-width:1100px; margin:0 auto; padding:24px 16px 64px; }}
  section {{ margin:26px 0; }}
  h2 {{ color:{T['deep']}; border-left:5px solid {T['matcha']}; padding-left:10px; }}
  .card {{ background:{T['panel']}; border:1px solid {T['light']}; border-radius:14px;
           padding:14px; box-shadow:0 2px 10px rgba(74,124,60,.08); }}
  .verdict {{ background:{T['panel']}; border:2px solid {T['matcha']}; border-radius:16px;
             padding:18px 22px; margin-top:18px; }}
  .verdict-title {{ color:{T['accent']}; font-weight:700; text-transform:uppercase;
                   letter-spacing:1px; font-size:13px; }}
  .verdict-main {{ font-size:24px; margin:6px 0; }}
  .verdict-sub {{ color:{T['deep']}; }}
  .verdict-basis {{ color:{T['stone']}; font-size:13px; margin-top:6px; }}
  .verdict-watch {{ margin-top:8px; color:#9a7d1f; font-size:13px; }}
  table.summary {{ border-collapse:collapse; width:100%; font-size:13px; }}
  table.summary th, table.summary td {{ border-bottom:1px solid {T['light']};
        padding:6px 10px; text-align:right; }}
  table.summary th {{ color:{T['deep']}; }}
  footer {{ text-align:center; color:{T['deep']}; padding:24px; font-size:13px; }}
</style></head>
<body>
  <header>
    <h1>MatchACell &mdash; cluster stability report</h1>
    <div class="scent">{MATCHA_SCENT}</div>
    {perf_strip}
  </header>
  <main>
    {verdict_html}
    <section><h2>Per-resolution summary</h2><div class="card">{summary_html}</div></section>
    {sections}
  </main>
  <footer>{MATCHA_SCENT} &middot; generated {datetime.now():%Y-%m-%d %H:%M}</footer>
</body></html>"""


# =========================================================================== #
# Helpers + CLI
# =========================================================================== #
def make_clustree(adata, leiden_keys, outdir, edge_weight_threshold=0.05):
    if not _HAS_CLUSTREE:
        log.warning("pyclustree not installed -> skipping clustree.")
        return
    log.info("Building clustree...")
    fig = clustree(adata, leiden_keys, title="Clustree (MatchACell)",
                   edge_weight_threshold=edge_weight_threshold, show_fraction=True)
    fig.set_size_inches(20, 20)
    fig.savefig(outdir / f"clustree_{edge_weight_threshold}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def sorted_leiden_keys(adata):
    keys = [c for c in adata.obs.columns if c.startswith("leiden_")]
    return sorted(keys, key=lambda k: float(k.split("_", 1)[1]))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="MatchACell step 1: data-driven QC + Leiden stability (GPU).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, type=Path, help="Input .h5ad (raw counts).")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    p.add_argument("--backend", choices=["auto", "gpu", "cpu"], default="auto")
    p.add_argument("--stability-backend", choices=["auto", "gpu", "cpu"], default="auto",
                   help="Backend for the bootstrap loop specifically. 'auto' runs "
                        "it on CPU for small datasets (see --gpu-cell-threshold), "
                        "where many tiny re-clusterings are faster on CPU.")
    p.add_argument("--gpu-cell-threshold", type=int, default=50000,
                   help="With --stability-backend auto, datasets smaller than this "
                        "(cells) run the bootstrap on CPU.")
    p.add_argument("--rmm-pool-fraction", type=float, default=0.5,
                   help="Fraction of free VRAM to pre-allocate for the RMM pool.")
    p.add_argument("--rmm-managed", action="store_true",
                   help="Use managed (unified) memory for the RMM pool.")
    p.add_argument("--skip-perf", action="store_true",
                   help="Disable performance profiling outputs.")
    p.add_argument("--no-color", action="store_true", help="Disable matcha CLI colours.")

    # clustering
    p.add_argument("--resolutions", type=float, nargs="+",
                   default=[0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    p.add_argument("--n-pcs", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--use-hvg", choices=["auto", "yes", "no"], default="auto")

    # embeddings
    p.add_argument("--skip-embeddings", action="store_true")
    p.add_argument("--skip-tsne", action="store_true", help="Skip t-SNE (the slowest).")
    p.add_argument("--skip-umap", action="store_true")

    # QC
    p.add_argument("--qc-nmads", type=float, default=5.0)
    p.add_argument("--qc-nmads-pct", type=float, default=3.0)
    p.add_argument("--min-counts", type=int, default=10)
    p.add_argument("--min-genes", type=int, default=5)
    p.add_argument("--mt-pct-hard", type=float, default=20.0)
    p.add_argument("--skip-qc", action="store_true")

    # stability
    p.add_argument("--n-iter", type=int, default=100)
    p.add_argument("--fraction", type=float, default=0.8)
    p.add_argument("--skip-stability", action="store_true")
    p.add_argument("--skip-transitions", action="store_true")
    p.add_argument("--skip-html", action="store_true")

    # verdict
    p.add_argument("--verdict-high", type=float, default=0.85)
    p.add_argument("--verdict-risk", type=float, default=0.60)
    p.add_argument("--verdict-min-median", type=float, default=0.75)
    p.add_argument("--verdict-max-risk-frac", type=float, default=0.25)

    p.add_argument("--edge-weight-threshold", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    Theme.configure(args.no_color)
    setup_logging()
    print_banner()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    rsc, use_gpu = init_backend(args.backend, pool_fraction=args.rmm_pool_fraction,
                                managed=args.rmm_managed)

    perf = None if args.skip_perf else PerfTracker(
        backend_label="gpu" if use_gpu else "cpu")

    def stage(name):
        return perf.stage(name) if perf is not None else contextlib.nullcontext()

    log.info("Loading %s", args.input)
    with stage("load"):
        adata = sc.read_h5ad(args.input)
    log.info("Loaded: %d cells x %d genes.", adata.n_obs, adata.n_vars)

    qc_obs = None
    if not args.skip_qc:
        with stage("qc"):
            adata, qc_summary, qc_obs = data_driven_qc(
                adata, nmads=args.qc_nmads, nmads_pct=args.qc_nmads_pct,
                min_counts_floor=args.min_counts, min_genes_floor=args.min_genes,
                mt_pct_hard=args.mt_pct_hard, outdir=args.outdir)
        (args.outdir / "qc").mkdir(parents=True, exist_ok=True)
        pd.Series(qc_summary).to_csv(args.outdir / "qc" / "qc_summary.csv")
    else:
        log.info("Skipping QC (--skip-qc).")

    use_hvg = (args.use_hvg == "yes") or (args.use_hvg == "auto" and adata.n_vars > 2000)
    log.info("HVG selection: %s", use_hvg)
    with stage("preprocess_cluster_embed"):
        adata = preprocess_and_cluster(
            adata, rsc, use_gpu, args.resolutions, args.n_pcs, args.n_neighbors,
            use_hvg, do_tsne=not args.skip_tsne, do_umap=not args.skip_umap)

    clustered_path = args.outdir / "clustered_multi_resolution.h5ad"
    adata.uns["MatchACell"] = {"scent": MATCHA_SCENT}
    adata.write(clustered_path)
    log.info("Saved clustered object -> %s", clustered_path)

    leiden_keys = sorted_leiden_keys(adata)
    with stage("clustree"):
        make_clustree(adata, leiden_keys, args.outdir, args.edge_weight_threshold)

    if not args.skip_embeddings:
        with stage("embedding_plots"):
            plot_embeddings(adata, leiden_keys, args.outdir)

    # Decide the bootstrap backend (it can differ from the main backend).
    stab_use_gpu = use_gpu
    if args.stability_backend == "cpu":
        stab_use_gpu = False
    elif args.stability_backend == "gpu":
        stab_use_gpu = use_gpu  # honoured only if a GPU is actually available
    else:  # auto
        if use_gpu and adata.n_obs < args.gpu_cell_threshold:
            log.info("Stability bootstrap -> CPU: %d cells < --gpu-cell-threshold "
                     "%d (many tiny re-clusterings run faster on CPU).",
                     adata.n_obs, args.gpu_cell_threshold)
            stab_use_gpu = False
    if stab_use_gpu and rsc is None:
        stab_use_gpu = False

    master = verdict = summary = transitions = None
    if not args.skip_stability:
        with stage("stability"):
            master = run_all_stability(adata, leiden_keys, rsc, stab_use_gpu, args.n_iter,
                                       args.fraction, args.n_neighbors, args.outdir, rng,
                                       perf=perf)
        verdict = compute_verdict(master, high=args.verdict_high, risk=args.verdict_risk,
                                  min_median=args.verdict_min_median,
                                  max_risk_frac=args.verdict_max_risk_frac)
        summary = verdict["summary"]
        (args.outdir / "MatchA_Verdict.txt").write_text(render_verdict_text(verdict),
                                                        encoding="utf-8")
        summary.round(4).to_csv(args.outdir / "stability" / "resolution_summary.csv")
        adata.uns["MatchACell"]["verdict_resolution"] = str(verdict["recommended"])
        adata.write(clustered_path)  # persist verdict into uns
        print_verdict(verdict)

    if not args.skip_transitions:
        with stage("transitions"):
            transitions = run_transitions(adata, leiden_keys, args.outdir)

    # Finalize performance metadata BEFORE assembling the HTML report so the
    # report's header strip carries the cell count and total runtime. (The
    # html_report stage itself is therefore not included in total_seconds; it
    # is still recorded as a stage and written by write_performance below.)
    if perf is not None:
        perf.finalize(adata.n_obs, adata.n_vars)

    if not args.skip_html:
        with stage("html_report"):
            build_html_report(adata, master, verdict, qc_obs, summary, transitions,
                              args.outdir, leiden_keys, perf=perf)

    if perf is not None:
        write_performance(perf, args.outdir)

    print(Theme.matcha(_TEALEAF))
    log.info("Done. Outputs in %s", args.outdir)
    print(Theme.deep("  " + MATCHA_SCENT))


if __name__ == "__main__":
    main()
