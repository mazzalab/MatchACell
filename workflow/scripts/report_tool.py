#!/usr/bin/env python3
"""
MatchACell - one HTML report for Steps 1-3 of a sample.

Reads the outputs of Step 1 (QC, Leiden sweep, bootstrap stability, Verdict),
Step 2 (the annotators' calls, as gathered by Step 3) and Step 3 (consensus)
and writes a single self-contained page: plotly.js is inlined, so the file
opens offline and can be shared as is.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as pyo
from plotly.subplots import make_subplots

SCRIPT_DIR = Path(__file__).resolve().parent

# Category colours: scanpy's default_20, so the report matches the PNG plots.
# Categories come most common first, so the most distinct colours go to them.
PALETTE_20 = ["#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2",
              "#b5bd61", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
              "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31"]
# One-hue sequential ramp (light -> dark) for magnitudes.
SEQ = [[0, "#e8f1fc"], [0.25, "#9ec5f4"], [0.5, "#3987e5"], [0.75, "#1c5cab"], [1, "#0d366b"]]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
CONSENSUS_STATUS = {"agreed": ("good", "✓"), "partial": ("warning", "◐"), "none": ("critical", "✕")}
RELATION = {"same": ("agree", "✓", "same term as the consensus"),
            "finer": ("agree", "✓", "a more specific term under the consensus"),
            "coarser": ("coarse", "≈", "a coarser term above the consensus"),
            "conflict": ("conflict", "✕", "a term outside the consensus lineage"),
            "no consensus": ("nocons", "–", "no consensus in this cluster")}
MAX_POINTS = 60_000
GRID = "rgba(128,128,128,0.22)"


def palette(n):
    """n category colours: scanpy's 20, then evenly spaced hues."""
    extra = [f"hsl({(i * 137.508) % 360:.0f},{60 if i % 2 else 45}%,{40 if i % 2 else 62}%)"
             for i in range(max(0, n - len(PALETTE_20)))]
    return (PALETTE_20 + extra)[:n]


def esc(x):
    return html.escape("" if x is None or (isinstance(x, float) and math.isnan(x)) else str(x))


def fmt(x, digits=2, pct=False):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    return f"{100 * x:.0f}%" if pct else f"{x:.{digits}f}"


def cluster_key(x):
    return [(0, int(p)) if p.isdigit() else (1, p) for p in str(x).split("_")]


def fig_html(fig, height=420):
    fig.update_layout(height=height, margin=dict(l=60, r=20, t=40, b=50),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12))
    fig.update_xaxes(gridcolor=GRID, zerolinecolor=GRID)
    fig.update_yaxes(gridcolor=GRID, zerolinecolor=GRID)
    return fig.to_html(full_html=False, include_plotlyjs=False,
                       config={"displaylogo": False, "responsive": True})


def table_html(df: pd.DataFrame, cls="data", raw_cols=(), num_cols=()):
    """Plain HTML table; `raw_cols` hold ready-made HTML, `num_cols` align right."""
    head = "".join(f'<th class="{"num" if c in num_cols else ""}">{esc(c)}</th>' for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c] if c in raw_cols else esc(r[c])
            cells.append(f'<td class="{"num" if c in num_cols else ""}">{v}</td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<div class="tablewrap"><table class="{cls}"><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def meter(x, kind=""):
    """Inline bar for a share in [0, 1] with its value."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    return (f'<span class="meter {kind}"><span style="width:{max(0, min(1, x)) * 100:.0f}%">'
            f'</span></span><span class="mv">{x:.2f}</span>')


def stability_class(x, risk=0.60, high=0.85):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return "risk" if x < risk else ("high" if x >= high else "")


def chip(status):
    kind, icon = CONSENSUS_STATUS.get(status, ("critical", "?"))
    return f'<span class="chip {kind}">{icon} {esc(status)}</span>'


def cl_link(cl_id, name):
    if isinstance(cl_id, str) and cl_id.startswith("CL:"):
        url = "http://purl.obolibrary.org/obo/" + cl_id.replace(":", "_")
        return f'<a href="{url}" target="_blank" rel="noopener">{esc(name)}</a> <span class="muted">{esc(cl_id)}</span>'
    return f'{esc(name)} <span class="muted">free text</span>'


def image_tag(path: Path, alt):
    import base64
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode()
    return f'<figure><img alt="{esc(alt)}" src="data:image/png;base64,{data}"><figcaption>{esc(alt)}</figcaption></figure>'


# =========================================================================== #
# Inputs
# =========================================================================== #
class Run:
    """Everything the report reads, with missing pieces left as None."""

    def __init__(self, step1: Path, consensus: Path, markers=None, settings=None):
        self.step1, self.cons_dir = step1, consensus
        self.settings = settings or {}
        a = ad.read_h5ad(step1 / "clustered_multi_resolution.h5ad", backed="r")
        self.obs = a.obs.copy()
        self.n_genes = a.n_vars
        self.basis = next((b for b in ("X_umap", "X_tsne", "X_pca") if b in a.obsm), None)
        self.xy = np.asarray(a.obsm[self.basis][:, :2]) if self.basis else None
        self.markers = self._marker_expression(a, markers) if markers else None
        a.file.close()
        self.leiden = sorted((c for c in self.obs if re.fullmatch(r"leiden_[\d.]+", c)),
                             key=lambda c: float(c.split("_")[1]))
        self.verdict = self._verdict(step1 / "MatchA_Verdict.txt")
        self.res_col = f"leiden_{self.verdict['resolution']}"
        self.qc = self._kv(step1 / "qc" / "qc_summary.csv")
        self.summary = self._csv(step1 / "stability" / "resolution_summary.csv",
                                 dtype={"Resolution": str})
        self.stability = {}
        for f in sorted((step1 / "stability" / "csv").glob("stability_*.csv")):
            df = pd.read_csv(f, dtype={"Cluster": str})
            self.stability[f.stem.removeprefix("stability_")] = df
        self.perf = self._csv(step1 / "performance" / "performance_metrics.csv")
        self.perf_meta = self._json(step1 / "performance" / "performance_summary.json")
        # Step 3 (and Step 2 through it)
        self.table = self._csv(consensus / "consensus_by_cluster.csv", dtype={"cluster": str})
        self.long = self._csv(consensus / "annotation_long.csv", dtype={"cluster": str})
        self.vocab = self._csv(consensus / "label_vocabulary.csv")
        self.pairs = self._csv(consensus / "annotator_agreement.csv")
        self.cons_settings = self._json(consensus / "consensus_settings.json") or {}
        self.validation = self._csv(consensus / "validation.csv")
        self.validation_by = self._csv(consensus / "validation_by_label.csv")
        calls = consensus / "cell_calls.csv.gz"
        self.calls = (pd.read_csv(calls, index_col="cell", dtype=str, keep_default_na=False)
                      if calls.exists() else None)

    @staticmethod
    def _csv(path, **kw):
        return pd.read_csv(path, **kw) if Path(path).exists() else None

    @staticmethod
    def _json(path):
        return json.loads(Path(path).read_text()) if Path(path).exists() else None

    @staticmethod
    def _kv(path):
        if not path.exists():
            return {}
        df = pd.read_csv(path, header=None, skiprows=1)
        return dict(zip(df[0], df[1]))

    @staticmethod
    def _verdict(path):
        text = path.read_text()

        def grab(pattern, default=None):
            m = re.search(pattern, text)
            return m.group(1).strip() if m else default
        watch = re.findall(r"cluster (\S+): Jaccard ([\d.]+)", text)
        return {
            "resolution": grab(r"Recommended resolution\s*:\s*leiden_(\S+)"),
            "n_clusters": grab(r"Number of clusters\s*:\s*(\d+)"),
            "median": grab(r"Median Jaccard\s*:\s*([\d.]+)"),
            "basis": grab(r"Selection basis\s*:\s*(.+)"),
            "high": float(grab(r"Clusters >= ([\d.]+)", 0.85)),
            "risk": float(grab(r"Clusters <\s+([\d.]+)", 0.60)),
            "edge": "largest resolution swept" in text,
            "watch": [(c, float(j)) for c, j in watch],
        }

    @staticmethod
    def _marker_expression(a, workbook):
        """Mean log-expression and fraction of expressing cells per cluster
        for the marker genes of the signature workbook."""
        try:
            sheets = pd.read_excel(workbook, sheet_name=None, header=None)
        except Exception as e:  # noqa: BLE001 - the report must not fail on it
            print(f"[report] marker workbook not read ({e}); skipping marker plot.")
            return None
        genes = {name: [g.strip() for g in df.iloc[:, 0].dropna().astype(str)
                        if g.strip() in a.var_names] for name, df in sheets.items()}
        genes = {k: list(dict.fromkeys(v)) for k, v in genes.items() if v}
        if not genes:
            return None
        flat = list(dict.fromkeys(g for v in genes.values() for g in v))
        idx = a.var_names.get_indexer(flat)
        X = a.layers["lognorm"] if "lognorm" in a.layers else a.X
        X = X[:, sorted(idx)] if hasattr(X, "shape") else X
        order = np.argsort(np.argsort(idx))
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        X = X[:, order]
        return {"genes": genes, "flat": flat, "X": X}


# =========================================================================== #
# Sections
# =========================================================================== #
def tiles(items):
    out = []
    for label, value, sub in items:
        out.append(f'<div class="tile"><div class="tl">{esc(label)}</div>'
                   f'<div class="tv">{value}</div><div class="ts">{sub}</div></div>')
    return f'<div class="tiles">{"".join(out)}</div>'


def overview(run: Run):
    v, t = run.verdict, run.table
    qc = run.qc
    items = []
    if qc:
        items.append(("Cells after QC", f"{int(qc['cells_out']):,}",
                      f"of {int(qc['cells_in']):,} ({qc['pct_removed']:.1f}% removed)"))
    else:
        items.append(("Cells", f"{len(run.obs):,}", "QC skipped (normalized input)"))
    items.append(("MatchA Verdict", f"leiden_{esc(v['resolution'])}",
                  f"{v['n_clusters']} clusters · median Jaccard {v['median']}"))
    if t is not None:
        settings = run.cons_settings.get("annotators", {})
        n_vot = sum(len(s["voters"]) for s in settings.values())
        items.append(("Annotators", f"{len(settings)}",
                      f"{n_vot} label sets · vote by {esc(run.cons_settings.get('vote_by', '?'))}"))
        counts = t["status"].value_counts()
        cells = t.groupby("status")["n_cells"].sum() / t["n_cells"].sum()
        items.append(("Consensus", f"{counts.get('agreed', 0)} / {len(t)} agreed",
                      " · ".join(f"{s}: {counts.get(s, 0)} clusters, {100 * cells.get(s, 0):.0f}% of cells"
                                 for s in ("partial", "none") if counts.get(s, 0))
                      or "every cluster agreed"))
    if run.validation is not None:
        c = run.validation.set_index("method").loc["Consensus"]
        items.append(("Against known labels", f"{100 * c['correct']:.0f}% correct",
                      f"{100 * c['coarser']:.0f}% coarser · {100 * c['wrong']:.0f}% wrong"))
    body = tiles(items)
    body += ('<p class="lead">The embedding below is shared by every step: colour it by a '
             'Leiden resolution, an annotator, the consensus or a QC metric. Click legend '
             'entries to hide or isolate them.</p>')
    body += '<div id="explorer" class="card"></div>'
    return body


def explorer_data(run: Run):
    """Data of the embedding explorer: coordinates once, then one colouring per option."""
    if run.xy is None:
        return None
    n = len(run.obs)
    rng = np.random.default_rng(0)
    keep = np.sort(rng.choice(n, MAX_POINTS, replace=False)) if n > MAX_POINTS else np.arange(n)
    cells = run.obs.index[keep]
    clusters = run.obs[run.res_col].astype(str).iloc[keep]
    options = []

    # colour of a CL term, shared by every annotator and the consensus
    term_of = {}
    if run.vocab is not None:
        for _, r in run.vocab.iterrows():
            term_of[(r["voter"], r["label"])] = r["cl_name"]
    if run.table is not None:
        for _, r in run.table.iterrows():
            term_of[("consensus", r["consensus_label"])] = r["consensus_cl_name"]
    term_weight = {}
    if run.calls is not None:
        for col in run.calls.columns.drop("cluster"):
            for lab, k in run.calls[col].value_counts().items():
                if lab:
                    t = term_of.get(("consensus" if col == "consensus_label" else col, lab), lab)
                    term_weight[t] = term_weight.get(t, 0) + k
    terms = sorted(term_weight, key=lambda t: -term_weight[t])
    term_color = dict(zip(terms, palette(len(terms))))

    def categorical(name, values, colors=None, group="Leiden"):
        values = pd.Series(values, index=cells).astype(str).replace({"nan": "", "None": ""})
        cats = [c for c in values.value_counts().index if c]
        if values.eq("").any():
            cats.append("")
        try:
            cats = sorted(cats, key=lambda c: (c == "", cluster_key(c))) if group == "Leiden" else cats
        except TypeError:
            pass
        code = {c: i for i, c in enumerate(cats)}
        cols = [colors(c) if colors else None for c in cats]
        if not colors:
            cols = palette(len([c for c in cats if c])) + (["#c8c8c8"] if "" in cats else [])
        cols = ["#c8c8c8" if c == "" else col for c, col in zip(cats, cols)]
        options.append({"name": name, "group": group, "kind": "cat",
                        "cats": [c or "no call" for c in cats], "colors": cols,
                        "codes": [code[v] for v in values]})

    def numeric(name, values, group, fmt_=".2f"):
        vals = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
        options.append({"name": name, "group": group, "kind": "num", "fmt": fmt_,
                        "values": [None if math.isnan(x) else round(float(x), 4) for x in vals]})

    if run.calls is not None:
        cc = run.calls.reindex(cells)
        categorical("Consensus", cc["consensus_label"].values,
                    lambda c: term_color.get(term_of.get(("consensus", c), c), "#c8c8c8"),
                    group="Step 3")
    for key in run.leiden:
        categorical(key + (" (Verdict)" if key == run.res_col else ""),
                    run.obs[key].astype(str).iloc[keep].values, group="Step 1 · Leiden")
    if run.calls is not None:
        for col in run.calls.columns.drop(["cluster", "consensus_label"]):
            categorical(col, cc[col].values,
                        lambda c, col=col: term_color.get(term_of.get((col, c), c), "#c8c8c8"),
                        group="Step 2 · annotators")
    if run.table is not None:
        t = run.table.set_index("cluster")
        numeric("Consensus agreement", clusters.map(t["agreement"]).values, "Step 3")
        if "stability" in t:
            numeric(f"Cluster stability (Jaccard, {run.res_col})",
                    clusters.map(t["stability"]).values, "Step 1 · stability")
    for col, label in [("total_counts", "Total counts"), ("n_genes_by_counts", "Detected genes"),
                       ("pct_counts_mt", "Mitochondrial %"), ("pct_counts_control", "Control probe %"),
                       ("doublet_score", "Doublet score (scDblFinder)")]:
        if col in run.obs:
            numeric(label, run.obs[col].iloc[keep].values, "Step 1 · QC",
                    ".2f" if col == "doublet_score" else (
                        ".0f" if "counts" in col and "pct" not in col or "genes" in col else ".1f"))

    xy = run.xy[keep]
    return {"x": np.round(xy[:, 0], 3).tolist(), "y": np.round(xy[:, 1], 3).tolist(),
            "cluster": clusters.tolist(), "resolution": run.res_col,
            "basis": run.basis.removeprefix("X_").upper(), "n_total": n, "n_shown": len(keep),
            "options": options, "seq": SEQ}


def step1(run: Run):
    v = run.verdict
    parts = []
    # --- QC ---------------------------------------------------------------
    qc = run.qc
    parts.append('<h3 id="s1-qc">Quality control</h3>')
    if qc:
        parts.append(
            '<p>Cells are removed when a metric is more than 5 median absolute deviations '
            '(MADs) from the median (counts, genes, top-20-gene share), more than 3 MADs '
            '(mitochondrial or control-probe %), above the hard mitochondrial ceiling, or '
            'below the absolute floors'
            + (', or called a doublet by scDblFinder' if "removed_by_doublet" in qc else '')
            + f'. <b>{int(qc["removed"]):,}</b> of '
            f'{int(qc["cells_in"]):,} cells ({qc["pct_removed"]:.1f}%) were removed; '
            f'{int(qc["cells_out"]):,} cells and {run.n_genes:,} genes go on.</p>')
        crit = {k.removeprefix("removed_by_"): v_ for k, v_ in qc.items()
                if k.startswith("removed_by_")}
        names = {"counts": "total counts", "genes": "detected genes", "mt": "mitochondrial %",
                 "control": "control probes", "floor": "absolute floor",
                 "doublet": "doublets (scDblFinder)"}
        f = go.Figure(go.Bar(x=[crit[k] for k in crit], y=[names.get(k, k) for k in crit],
                             orientation="h", marker_color="#3987e5",
                             hovertemplate="%{y}: %{x:,} cells<extra></extra>",
                             text=[f"{int(crit[k]):,}" for k in crit], textposition="outside"))
        f.update_layout(xaxis_title="cells flagged (a cell can fail several criteria)",
                        yaxis=dict(autorange="reversed"), showlegend=False)
        parts.append(f'<div class="card">{fig_html(f, 260)}</div>')
        imgs = "".join(image_tag(run.step1 / "qc" / name, cap) for name, cap in [
            ("qc_hist_mad_bounds.png", "QC metrics with the MAD bounds"),
            ("qc_violin_kept_vs_removed.png", "Kept vs removed cells"),
            ("qc_counts_vs_genes.png", "Counts vs detected genes"),
            ("qc_barcode_rank.png", "Barcode rank (knee) plot"),
            ("qc_highest_expr_genes.png", "Most expressed genes")])
        if imgs:
            parts.append(f'<details><summary>QC diagnostic figures</summary><div class="figs">{imgs}</div></details>')
    else:
        parts.append('<p>The input was already normalized (<code>type: normalized</code>): '
                     'QC and preprocessing were skipped.</p>')

    # --- Verdict ------------------------------------------------------------
    parts.append('<h3 id="s1-verdict">Resolution sweep and the MatchA Verdict</h3>')
    parts.append(
        f'<p>Leiden clustering was run at {len(run.leiden)} resolutions. Each cluster\'s '
        f'stability is its mean best-match Jaccard index over bootstrap subsamples of 80% '
        f'of the cells (Hennig 2007): ≥ {v["high"]:.2f} is highly stable, &lt; '
        f'{v["risk"]:.2f} is at risk. The Verdict picks the finest resolution whose median '
        f'stability is ≥ 0.75, with at most 25% of its clusters and 5% of its cells in '
        f'clusters at risk.</p>')
    watch = ", ".join(f"{esc(c)} ({j:.2f})" for c, j in v["watch"]) or "none"
    edge = ('<div class="warnline">⚠ This is the largest resolution swept: the data may stay '
            'stable at higher resolutions. Extend <code>--resolutions</code> or check a coarser '
            'stable resolution.</div>' if v["edge"] else "")
    parts.append(f'''<div class="verdict"><div class="vt">MatchA Verdict</div>
      <div class="vm">Annotate at <b>leiden_{esc(v["resolution"])}</b></div>
      <div class="vs">{esc(v["n_clusters"])} clusters · median Jaccard {esc(v["median"])} ·
      {esc(v["basis"])}</div><div class="vw">Clusters at risk: {watch}</div>{edge}</div>''')
    if run.summary is not None:
        s = run.summary.copy()
        rec = v["resolution"]
        s["Resolution"] = s["Resolution"].astype(str)
        disp = pd.DataFrame({
            "resolution": [f"<b>{esc(r)} ◀ Verdict</b>" if r == rec else esc(r) for r in s["Resolution"]],
            "clusters": s["n_clusters"],
            "median Jaccard": s["median"].map(lambda x: fmt(x)),
            "mean": s["mean"].map(lambda x: fmt(x)),
            "min": s["min"].map(lambda x: fmt(x)),
            f"≥ {v['high']:.2f}": s["frac_high"].map(lambda x: fmt(x, pct=True)),
            f"< {v['risk']:.2f}": s["frac_risk"].map(lambda x: fmt(x, pct=True)),
        })
        if "frac_risk_cells" in s:
            disp[f"cells in < {v['risk']:.2f}"] = s["frac_risk_cells"].map(
                lambda x: f"{100 * x:.1f}%")
        parts.append(table_html(disp, raw_cols={"resolution"},
                                num_cols=set(disp.columns) - {"resolution"}))
    if run.stability:
        order = sorted(run.stability, key=float)
        f = go.Figure()
        for r in order:
            df = run.stability[r]
            is_rec = r == v["resolution"]
            f.add_trace(go.Box(
                y=df["Jaccard_Stability"], name=r, boxpoints="all", jitter=0.5, pointpos=0,
                marker=dict(size=6, color="#1c5cab" if is_rec else "#86b6ef",
                            line=dict(width=1, color="white")),
                line=dict(color="#1c5cab" if is_rec else "#9a9890", width=1.5),
                fillcolor="rgba(0,0,0,0)", text=[f"cluster {c}" for c in df["Cluster"]],
                hovertemplate="resolution " + r + "<br>%{text}: %{y:.3f}<extra></extra>"))
        f.add_hline(y=v["high"], line_dash="dot", line_color="#0ca30c",
                    annotation_text=f"high {v['high']:.2f}", annotation_position="top left")
        f.add_hline(y=v["risk"], line_dash="dot", line_color="#d03b3b",
                    annotation_text=f"risk {v['risk']:.2f}", annotation_position="bottom left")
        f.update_layout(showlegend=False, xaxis_title="Leiden resolution (dark: the Verdict)",
                        yaxis_title="cluster stability (Jaccard)", yaxis_range=[0, 1.04])
        parts.append(f'<div class="card">{fig_html(f, 400)}</div>')

    # --- per-cluster stability at the Verdict --------------------------------
    df = run.stability.get(v["resolution"])
    if df is not None:
        df = df.assign(_k=df["Cluster"].map(cluster_key)).sort_values("_k")
        colors = ["#d03b3b" if x < v["risk"] else ("#1c5cab" if x >= v["high"] else "#86b6ef")
                  for x in df["Jaccard_Stability"]]
        sizes = run.obs[run.res_col].astype(str).value_counts()
        f = go.Figure(go.Bar(x=df["Cluster"], y=df["Jaccard_Stability"], marker_color=colors,
                             customdata=df["Cluster"].map(sizes).fillna(0).astype(int),
                             hovertemplate="cluster %{x}: %{y:.3f}<br>%{customdata:,} cells<extra></extra>"))
        f.add_hline(y=v["risk"], line_dash="dot", line_color="#d03b3b")
        f.update_layout(xaxis_title=f"cluster ({run.res_col}); red: below the risk threshold",
                        yaxis_title="stability (Jaccard)", yaxis_range=[0, 1.04],
                        xaxis_type="category")
        parts.append(f'<div class="card">{fig_html(f, 300)}</div>')

    # --- cluster flow across resolutions -----------------------------------
    if len(run.leiden) > 1:
        parts.append('<h3 id="s1-flow">How clusters split across resolutions</h3>'
                     '<p>Each column is a resolution, each block a cluster; bands show where '
                     'the cells of a cluster go at the next resolution (bands under 5% of a '
                     'cluster are hidden). Dark blocks: the Verdict resolution.</p>')
        parts.append(f'<div class="card">{fig_html(sankey(run), 520)}</div>')

    # --- performance ----------------------------------------------------------
    if run.perf is not None:
        p = run.perf
        meta = run.perf_meta or {}
        f = go.Figure(go.Bar(x=p["seconds"], y=p["stage"], orientation="h", marker_color="#86b6ef",
                             hovertemplate="%{y}: %{x:.1f} s<extra></extra>"))
        f.update_layout(xaxis_title="seconds", yaxis=dict(autorange="reversed"))
        parts.append(f'<details><summary>Run time of Step 1 ({meta.get("backend", "?")} backend, '
                     f'total {meta.get("total_seconds", float("nan")):.0f} s)</summary>'
                     f'<div class="card">{fig_html(f, 280)}</div></details>')
    return "\n".join(parts)


def sankey(run: Run):
    cols = run.leiden
    labels, index, colors, xs = [], {}, [], []
    for i, col in enumerate(cols):
        res = col.split("_")[1]
        for c in sorted(run.obs[col].astype(str).unique(), key=cluster_key):
            index[(col, c)] = len(labels)
            labels.append(f"{res}: {c}")
            colors.append("#1c5cab" if col == run.res_col else "#b7d3f6")
            xs.append(0.001 + 0.998 * i / max(1, len(cols) - 1))
    src, dst, val = [], [], []
    for a, b in zip(cols[:-1], cols[1:]):
        ct = pd.crosstab(run.obs[a].astype(str), run.obs[b].astype(str))
        for pa in ct.index:
            for pb in ct.columns:
                n = int(ct.at[pa, pb])
                if n and n >= 0.05 * ct.loc[pa].sum():
                    src.append(index[(a, pa)])
                    dst.append(index[(b, pb)])
                    val.append(n)
    f = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(label=labels, color=colors, pad=6, thickness=12, x=xs,
                  line=dict(width=0),
                  hovertemplate="resolution %{label}<br>%{value:,} cells<extra></extra>"),
        link=dict(source=src, target=dst, value=val, color="rgba(57,135,229,0.18)",
                  hovertemplate="%{source.label} → %{target.label}: %{value:,} cells<extra></extra>")))
    return f


def calls_grid(run: Run):
    """Clusters x voters: each voter's leading label, coloured by how it relates
    to the cluster's consensus."""
    long, t = run.long, run.table
    voters = list(dict.fromkeys(long["voter"]))
    lead = (long.sort_values("share", ascending=False)
                .drop_duplicates(["cluster", "voter"]).set_index(["cluster", "voter"]))
    head = "".join(f"<th>{esc(v)}</th>" for v in voters)
    rows = []
    for _, r in t.iterrows():
        cl = r["cluster"]
        cells = [f'<td class="num">{esc(cl)}</td>',
                 f'<td>{chip(r["status"])} {esc(r["consensus_label"])}</td>']
        for v in voters:
            if (cl, v) not in lead.index:
                cells.append('<td class="rel nocall">no call</td>')
                continue
            x = lead.loc[(cl, v)]
            kind, icon, why = RELATION.get(x["relation"], ("nocons", "?", ""))
            cells.append(f'<td class="rel {kind}" title="{esc(x["cl_name"])}: {esc(why)}">'
                         f'{icon} {esc(x["label"])} <span class="muted">{100 * x["share"]:.0f}%</span></td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    legend = " ".join(f'<span class="rel {k} key">{i} {esc(d)}</span>' for k, i, d in [
        ("agree", "✓", "same term or a more specific one"),
        ("coarse", "≈", "coarser, compatible"), ("conflict", "✕", "conflicting lineage"),
        ("nocall", "", "no call (Unknown)")])
    return (f'<div class="legend">{legend}</div><div class="tablewrap"><table class="data grid">'
            f'<thead><tr><th>cluster</th><th>consensus</th>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def step2(run: Run):
    if run.long is None:
        return "<p>Step 3 outputs are missing: run the <code>consensus</code> target.</p>"
    parts = []
    st = run.cons_settings.get("annotators", {})
    vocab = run.vocab
    rows = []
    n = len(run.obs)
    for a, s in st.items():
        for v in s["voters"]:
            vv = vocab[vocab["voter"] == v]
            called = vv["n_cells"].sum()
            how = vv.groupby("mapped_by")["n_cells"].sum() / max(called, 1)
            rows.append({
                "annotator": a + (f" · {v.split(':', 1)[1]}" if ":" in v else ""),
                "family": s["family"],
                "labels": len(vv),
                "cells with a call": f"{100 * called / n:.0f}%",
                "mapped to CL": f"{100 * (1 - how.get('unmapped', 0)):.0f}%",
                "weight": fmt(s["weight"], 2),
            })
    parts.append('<h3 id="s2-annotators">Annotators</h3>')
    parts.append('<p>Each annotator labels the cells on its own. Marker-signature methods score '
                 'the gene lists of the signature workbook; reference models transfer labels '
                 'from annotated atlases; LLM annotators name clusters from their marker genes. '
                 'Every label is then mapped to a Cell Ontology (CL) term for the vote.</p>')
    d = pd.DataFrame(rows)
    parts.append(table_html(d, num_cols={"labels", "cells with a call", "mapped to CL", "weight"}))
    parts.append('<h3 id="s2-calls">Calls per cluster</h3>'
                 '<p>The leading label of every annotator in every cluster, with the share of '
                 'the cluster\'s cells that carry it. Colour and symbol show how the label '
                 'relates to the cluster\'s consensus (Step 3) on the Cell Ontology. Hover a '
                 'cell for the CL term.</p>')
    parts.append(calls_grid(run))
    if run.pairs is not None and len(run.pairs):
        parts.append('<h3 id="s2-agreement">How often the annotators agree</h3>'
                     '<p>Over the cells two annotators both label: the share with the same CL '
                     'term (left), and with compatible terms, where one is the other or an '
                     'ancestor of it, like "T cell" and "CD4-positive, alpha-beta T cell" '
                     '(right). Annotators that agree far more with each other than with the rest '
                     'share their evidence; that is why Step 3 can pool them by family.</p>')
        parts.append(f'<div class="card">{fig_html(agreement_heatmaps(run.pairs), 430)}</div>')
    sig = signature_heatmaps(run)
    if sig:
        parts.append(sig)
    if vocab is not None:
        v = vocab.copy()
        v["n_cells"] = v["n_cells"].fillna(0).astype(int)
        v = v.sort_values(["mapped_by", "voter", "n_cells"], ascending=[False, True, False])
        un = (v["mapped_by"] == "unmapped").sum()
        v["CL term"] = [cl_link(i, nme) for i, nme in zip(v["cl_id"], v["cl_name"])]
        v = v[["voter", "label", "CL term", "mapped_by", "n_cells"]].rename(
            columns={"voter": "annotator", "mapped_by": "mapped by", "n_cells": "cells"})
        flag = (f' · <b class="bad">{un} unmapped</b>: add them to a label map'
                if un else " · all mapped")
        parts.append(f'<details><summary>Label vocabulary: {len(v)} labels{flag}</summary>'
                     f'{table_html(v, raw_cols={"CL term"}, num_cols={"cells"})}</details>')
    return "\n".join(parts)


def agreement_heatmaps(pairs):
    voters = list(dict.fromkeys(list(pairs["voter_a"]) + list(pairs["voter_b"])))
    fig = make_subplots(rows=1, cols=2, subplot_titles=("same CL term", "compatible terms"),
                        horizontal_spacing=0.18)
    for k, col in enumerate(["same", "compatible"], start=1):
        m = pd.DataFrame(np.eye(len(voters)), index=voters, columns=voters)
        for _, r in pairs.iterrows():
            m.at[r["voter_a"], r["voter_b"]] = m.at[r["voter_b"], r["voter_a"]] = r[col]
        fig.add_trace(go.Heatmap(z=m.values, x=voters, y=voters, zmin=0, zmax=1,
                                 colorscale=SEQ, showscale=(k == 2),
                                 text=[[f"{100 * x:.0f}%" for x in row] for row in m.values],
                                 texttemplate="%{text}",
                                 hovertemplate="%{y} vs %{x}: %{text}<extra></extra>"),
                      row=1, col=k)
        fig.update_yaxes(autorange="reversed", row=1, col=k)
    return fig


def signature_heatmaps(run: Run):
    out = []
    for a, s in run.cons_settings.get("annotators", {}).items():
        if a not in ("ScoreGenes", "AddModuleScore"):
            continue
        xlsx = Path(s["file"]).parent / "cluster_annotation_summary.xlsx"
        if not xlsx.exists():
            continue
        df = pd.read_excel(xlsx, dtype={"cluster": str}).set_index("cluster")
        df = df.loc[sorted(df.index, key=cluster_key)]
        sigs = [c for c in df.columns if c not in ("final_annotation", "highest_score")]
        z = df[sigs].astype(float)
        lim = float(np.nanmax(np.abs(z.values))) or 1.0
        f = go.Figure(go.Heatmap(
            z=z.values, x=sigs, y=[f"{c}" for c in z.index], zmid=0, zmin=-lim, zmax=lim,
            colorscale=[[0, "#1c5cab"], [0.5, "#f0efec"], [1, "#d03b3b"]],
            text=[[f"{x:.2f}" for x in row] for row in z.values], texttemplate="%{text}",
            customdata=[[df.at[c, "final_annotation"]] * len(sigs) for c in z.index],
            hovertemplate="cluster %{y} · %{x}: %{z:.3f}<br>call: %{customdata}<extra></extra>"))
        f.update_layout(xaxis_title="signature", yaxis_title="cluster",
                        yaxis=dict(autorange="reversed", type="category"))
        out.append(f'<details><summary>{esc(a)}: mean signature score per cluster</summary>'
                   f'<p class="muted">The call is the highest-scoring signature above the '
                   f'threshold; the other signatures above it are listed after it in the '
                   f'annotator\'s own table.</p><div class="card">{fig_html(f, 30 * len(z) + 140)}</div></details>')
    return "".join(out)


def step3(run: Run):
    t = run.table
    if t is None:
        return "<p>No consensus yet: run the <code>consensus</code> target.</p>"
    s = run.cons_settings
    v = run.verdict
    parts = []
    parts.append(
        f'<p>Every annotator votes in each cluster with the share of the cluster\'s cells '
        f'carrying each label. A vote for a CL term also supports all its ancestors, so '
        f'"CD4-positive, alpha-beta T cell" and "CD8-positive, alpha-beta T cell" disagree '
        f'but both support "T cell" (the ontology vote of popV, Ergen et al. 2024). The '
        f'consensus is the most specific term supported by more than '
        f'<b>{100 * s.get("min_agreement", 0.5):.0f}%</b> of the vote. '
        + ('Each family of annotators (marker signatures, reference models, LLMs) has one '
           'vote, split between its members that labelled the cluster. '
           if s.get("vote_by") == "family" else 'Each annotator has one vote. ')
        + '<b>agreed</b>: some annotator gave that very term; <b>partial</b>: they only '
        'agree on a coarser common ancestor; <b>none</b>: no term reaches the bar.</p>')
    rows = []
    for _, r in t.iterrows():
        rows.append({
            "cluster": esc(r["cluster"]),
            "cells": f'{int(r["n_cells"]):,}',
            "stability": meter(r.get("stability"), stability_class(r.get("stability"), v["risk"], v["high"])),
            "consensus": f'<b>{esc(r["consensus_label"])}</b>',
            "Cell Ontology term": cl_link(r["consensus_cl_id"], r["consensus_cl_name"]),
            "agreement": meter(r["agreement"]),
            "compatible": meter(r["compatibility"]),
            "status": chip(r["status"]),
            "most voted single term": (esc(r["top_label"]) if r["status"] != "agreed"
                                       and r["top_label"] != r["consensus_cl_name"] else ""),
        })
    d = pd.DataFrame(rows)
    parts.append('<h3 id="s3-table">Consensus per cluster</h3>')
    parts.append('<p class="muted"><b>agreement</b>: share of the vote that supports the '
                 'consensus term or a finer one. <b>compatible</b>: adds the votes for coarser '
                 'terms above it, which do not contradict it. <b>stability</b>: the Step 1 '
                 'bootstrap Jaccard of the cluster (red below the risk threshold).</p>')
    parts.append(table_html(d, raw_cols=set(d.columns), num_cols={"cells"}))

    parts.append('<h3 id="s3-trust">Which clusters to trust</h3>'
                 '<p>Stability (Step 1) against agreement (Step 3), one point per cluster sized '
                 'by its cells. Top right: stable clusters the annotators agree on. Low '
                 'stability: the cluster may not be a real population at this resolution. '
                 'Low agreement on a stable cluster: a real population the annotators '
                 'disagree on, so check its markers below.</p>')
    parts.append(f'<div class="card">{fig_html(trust_map(run), 440)}</div>')

    dot = marker_dotplot(run)
    if dot:
        parts.append('<h3 id="s3-markers">Marker evidence</h3>'
                     '<p>Marker genes of the signature workbook in each cluster: dot size is the '
                     'share of cells expressing the gene, colour its mean log-normalized '
                     'expression scaled per gene. Rows are labelled with the consensus.</p>')
        parts.append(dot)
    if run.validation is not None:
        parts.append(validation_section(run))
    return "\n".join(parts)


def trust_map(run: Run):
    t = run.table.copy()
    v = run.verdict
    f = go.Figure()
    t["stability"] = pd.to_numeric(t.get("stability"), errors="coerce")
    for status, (kind, icon) in CONSENSUS_STATUS.items():
        d = t[t["status"] == status]
        if d.empty:
            continue
        f.add_trace(go.Scatter(
            x=d["stability"], y=d["agreement"], mode="markers+text", name=f"{icon} {status}",
            text=d["cluster"], textposition="top center",
            marker=dict(size=np.sqrt(d["n_cells"]) / np.sqrt(t["n_cells"].max()) * 40 + 8,
                        color=STATUS[kind], opacity=0.85, line=dict(width=2, color="white"),
                        symbol={"agreed": "circle", "partial": "diamond", "none": "x"}[status]),
            customdata=np.stack([d["consensus_label"], d["n_cells"]], axis=1),
            hovertemplate=("cluster %{text}: %{customdata[0]}<br>stability %{x:.2f} · "
                           "agreement %{y:.2f}<br>%{customdata[1]:,} cells<extra></extra>")))
    f.add_vline(x=v["risk"], line_dash="dot", line_color="#9a9890")
    f.add_hline(y=run.cons_settings.get("min_agreement", 0.5), line_dash="dot", line_color="#9a9890")
    f.update_layout(xaxis_title="Step 1 stability (Jaccard)", yaxis_title="Step 3 agreement",
                    xaxis_range=[0, 1.05], yaxis_range=[-0.02, 1.15],
                    legend=dict(orientation="h", y=1.1))
    return f


def marker_dotplot(run: Run):
    m = run.markers
    if not m or run.table is None:
        return ""
    clusters = run.obs[run.res_col].astype(str).values
    order = list(run.table["cluster"])
    cons = run.table.set_index("cluster")["consensus_label"]
    X = m["X"]
    mean = pd.DataFrame({c: X[clusters == c].mean(axis=0) for c in order}, index=m["flat"]).T
    frac = pd.DataFrame({c: (X[clusters == c] > 0).mean(axis=0) for c in order}, index=m["flat"]).T
    # up to 6 genes per signature, the most variable across clusters
    genes, groups = [], []
    for sig, gl in m["genes"].items():
        top = mean[gl].var().sort_values(ascending=False).index[:6]
        for g in top:
            if g not in genes:
                genes.append(g)
                groups.append(sig)
    scaled = (mean[genes] - mean[genes].min()) / (mean[genes].max() - mean[genes].min()).replace(0, 1)
    xs, ys, size, color, hover = [], [], [], [], []
    ylab = [f"{c} · {cons[c]}" for c in order]
    for i, c in enumerate(order):
        for j, g in enumerate(genes):
            xs.append(j)
            ys.append(ylab[i])
            size.append(float(frac.at[c, g]))
            color.append(float(scaled.at[c, g]))
            hover.append(f"{g} ({groups[j]}) in cluster {c}<br>{100 * frac.at[c, g]:.0f}% of cells, "
                         f"mean {mean.at[c, g]:.2f}")
    f = go.Figure(go.Scatter(
        x=xs, y=ys, mode="markers", text=hover, hovertemplate="%{text}<extra></extra>",
        marker=dict(size=[3 + 15 * s for s in size], color=color, colorscale=SEQ, cmin=0, cmax=1,
                    line=dict(width=0.5, color="rgba(128,128,128,0.5)"),
                    colorbar=dict(title="scaled<br>mean", thickness=10))))
    ticks = [f"{g}" for g in genes]
    f.update_xaxes(tickvals=list(range(len(genes))), ticktext=ticks, tickangle=-60)
    f.update_yaxes(autorange="reversed", type="category")
    # signature names above their genes
    start = 0
    for sig in dict.fromkeys(groups):
        k = groups.count(sig)
        f.add_annotation(x=start + (k - 1) / 2, y=1.0, yref="paper", yanchor="bottom",
                         text=sig, showarrow=False, font=dict(size=11))
        if start:
            f.add_vline(x=start - 0.5, line_color=GRID)
        start += k
    f.update_layout(margin=dict(t=60))
    return f'<div class="card">{fig_html(f, 26 * len(order) + 200)}</div>'


def validation_section(run: Run):
    val = run.validation.copy()
    col = run.cons_settings.get("truth_col")
    parts = [f'<h3 id="s3-validation">Against known labels (<code>{esc(col)}</code>)</h3>',
             '<p>Each cell\'s label from an annotator, or from its cluster\'s consensus, compared '
             'with the known label on the Cell Ontology. <b>correct</b>: the same term or a more '
             'specific one; <b>coarser</b>: an ancestor of the known term, right but less '
             'specific; <b>wrong</b>: another lineage. <b>Cluster ceiling</b> labels each cluster '
             'with its most common known label: no cluster-level annotation can do better at '
             'this resolution. Known labels are only as good as the way they were made: labels '
             'from marker genes favour marker methods.</p>']
    f = go.Figure()
    val["correct"] = val["exact"] + val["finer"]
    for key, name, color in [("correct", "correct", STATUS["good"]),
                             ("coarser", "coarser", STATUS["warning"]),
                             ("wrong", "wrong", STATUS["critical"]),
                             ("no call", "no call", "#c3c2b7")]:
        f.add_trace(go.Bar(y=val["method"], x=val[key], name=name, orientation="h",
                           marker=dict(color=color, line=dict(width=1, color="white")),
                           text=[f"{100 * x:.0f}%" if x >= 0.04 else "" for x in val[key]],
                           textposition="inside", insidetextanchor="middle",
                           hovertemplate="%{y}: %{x:.1%} " + name + "<extra></extra>"))
    f.update_layout(barmode="stack", xaxis=dict(tickformat=".0%", range=[0, 1]),
                    yaxis=dict(autorange="reversed"),
                    legend=dict(orientation="h", y=1.12, traceorder="normal"))
    parts.append(f'<div class="card">{fig_html(f, 90 + 42 * len(val))}</div>')
    if run.validation_by is not None:
        b = run.validation_by
        m = b.pivot(index="known_label", columns="method", values="correct")
        m = m[[c for c in val["method"] if c in m.columns]]
        n = b.drop_duplicates("known_label").set_index("known_label")["n_cells"]
        m = m.loc[n.sort_values(ascending=False).index]
        f = go.Figure(go.Heatmap(z=m.values, x=list(m.columns),
                                 y=[f"{i} ({n[i]:,})" for i in m.index], zmin=0, zmax=1,
                                 colorscale=SEQ, text=[[f"{100 * x:.0f}%" for x in r] for r in m.values],
                                 texttemplate="%{text}",
                                 hovertemplate="%{y} · %{x}: %{z:.1%} correct<extra></extra>"))
        f.update_layout(yaxis=dict(autorange="reversed"), xaxis_side="top")
        parts.append('<p class="muted">Share correct per known label (cells in brackets).</p>')
        parts.append(f'<div class="card">{fig_html(f, 30 * len(m) + 120)}</div>')
    return "\n".join(parts)


def methods(run: Run, settings):
    cs = run.cons_settings
    rows = [("Step 1", "QC", "MAD outliers: 5 MADs on log counts, log genes and the top-20-gene "
             "share; 3 MADs and a hard ceiling on mitochondrial %; Xenium control probes dropped."),
            ("Step 1", "Clustering", "normalize_total → log1p → HVG (PCA only) → PCA → kNN → Leiden "
             f"at {', '.join(k.split('_')[1] for k in run.leiden)}; embedding {esc(run.basis or '–')}."),
            ("Step 1", "Stability", "Bootstrap of 80% subsamples, Leiden re-run on the subsample's "
             "kNN graph, best-match Jaccard per cluster (Hennig 2007)."),
            ("Step 2", "Annotators", ", ".join(f"{esc(a)} ({esc(s['family'])})"
                                                for a, s in cs.get("annotators", {}).items()) or "–"),
            ("Step 3", "Ontology", f"Cell Ontology {esc(cs.get('ontology_version', '?'))}, "
             f"{cs.get('n_terms', '?')} terms, is_a only"),
            ("Step 3", "Label maps", ", ".join(esc(Path(p).name) for p in cs.get("label_maps", [])) or "none"),
            ("Step 3", "Vote", f"by {esc(cs.get('vote_by', '?'))}, more than "
             f"{100 * cs.get('min_agreement', 0.5):.0f}% needed")]
    body = table_html(pd.DataFrame(rows, columns=["step", "what", "how"]))
    if settings:
        body += ('<details><summary>Configuration of this run</summary><pre>'
                 + esc(json.dumps(settings, indent=2)) + "</pre></details>")
    refs = [
        "Germain P-L, Sonrel A, Robinson MD (2020). pipeComp, a general framework for the evaluation of computational pipelines, reveals performant single cell RNA-seq preprocessing tools. <i>Genome Biology</i> 21:227.",
        "Heumos L, Schaar AC, et al. (2023). Best practices for single-cell analysis across modalities. <i>Nature Reviews Genetics</i> 24:550–572.",
        "Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing well-connected communities. <i>Scientific Reports</i> 9:5233.",
        "Hennig C (2007). Cluster-wise assessment of cluster stability. <i>Computational Statistics &amp; Data Analysis</i> 52:258–271.",
        "Tirosh I, et al. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. <i>Science</i> 352:189–196 (gene-signature scores).",
        "Domínguez Conde C, et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. <i>Science</i> 376:eabl5197 (CellTypist).",
        "Xu C, et al. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. <i>Molecular Systems Biology</i> 17:e9620 (scANVI).",
        "Ergen C, et al. (2024). Consensus prediction of cell type labels in single-cell data with popV. <i>Nature Genetics</i> 56:2731–2738.",
        "Diehl AD, et al. (2016). The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability. <i>Journal of Biomedical Semantics</i> 7:44.",
        "Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene expression data analysis. <i>Genome Biology</i> 19:15.",
    ]
    body += "<h3>References</h3><ol class='refs'>" + "".join(f"<li>{r}</li>" for r in refs) + "</ol>"
    return body


def files_section(run: Run):
    rows = [
        ("clustered_multi_resolution.h5ad / .rds", "Step 1", "QC-filtered cells, every Leiden resolution, embeddings"),
        ("MatchA_Verdict.txt", "Step 1", "the recommended resolution and the clusters at risk"),
        ("MatchACell_report.html", "Step 1", "the Step 1 report"),
        ("qc/, embeddings/, stability/, transitions/, performance/", "Step 1", "figures and tables"),
        ("annotation/<annotator>/", "Step 2", "each annotator's labelled object, tables and figures"),
        ("consensus/consensus_by_cluster.csv", "Step 3", "the consensus of every cluster (the table above)"),
        ("consensus/consensus_annotated.h5ad", "Step 3", "the Step 1 object with the consensus in .obs"),
        ("consensus/annotation_long.csv", "Step 3", "every label of every annotator per cluster, its CL term and relation to the consensus"),
        ("consensus/label_vocabulary.csv", "Step 3", "every label, its CL term and how it was mapped"),
        ("consensus/annotator_agreement.csv", "Step 3", "pairwise agreement of the annotators"),
        ("consensus/cell_calls.csv.gz", "Step 3", "per-cell labels of every annotator and the consensus"),
        ("consensus/validation*.csv", "Step 3", "scores against known labels, when truth_col is set"),
    ]
    return table_html(pd.DataFrame(rows, columns=["file", "step", "content"]))


# =========================================================================== #
# Page
# =========================================================================== #
CSS = """
:root{--bg:#f6f7f2;--panel:#ffffff;--ink:#1f2a17;--ink2:#4c5a40;--muted:#7d8873;--line:#dfe6d3;
 --accent:#4a7c3c;--accent2:#7ab661;--good:#0ca30c;--warn:#b07a00;--bad:#c22f2f;
 --agree-bg:#e3f4e1;--coarse-bg:#fdf3d7;--conflict-bg:#fbe3e1;--nocall-bg:#f0efec;--meter:#3987e5}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#121410;--panel:#1b1e18;
 --ink:#eef1e8;--ink2:#c3c9b8;--muted:#8d977f;--line:#2f352a;--accent:#8fc77a;--accent2:#7ab661;
 --good:#3fc43f;--warn:#e0a920;--bad:#ef6b6b;--agree-bg:#1d3a1d;--coarse-bg:#3a3214;
 --conflict-bg:#432020;--nocall-bg:#262923;--meter:#5598e7}}
:root[data-theme="dark"]{--bg:#121410;--panel:#1b1e18;--ink:#eef1e8;--ink2:#c3c9b8;--muted:#8d977f;
 --line:#2f352a;--accent:#8fc77a;--accent2:#7ab661;--good:#3fc43f;--warn:#e0a920;--bad:#ef6b6b;
 --agree-bg:#1d3a1d;--coarse-bg:#3a3214;--conflict-bg:#432020;--nocall-bg:#262923;--meter:#5598e7}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{background:linear-gradient(135deg,#7ab661,#4a7c3c);color:#fff;padding:22px 24px}
header h1{margin:0;font-size:22px;font-weight:650}
header .sub{opacity:.92;font-size:13px;margin-top:4px}
nav{position:sticky;top:0;z-index:5;background:var(--panel);border-bottom:1px solid var(--line);
 display:flex;gap:4px;padding:6px 16px;overflow-x:auto}
nav a{color:var(--ink2);text-decoration:none;padding:6px 10px;border-radius:8px;white-space:nowrap;font-size:13px}
nav a:hover{background:var(--nocall-bg);color:var(--ink)}
nav .sp{flex:1}
nav button{background:none;border:1px solid var(--line);color:var(--ink2);border-radius:8px;padding:4px 10px;cursor:pointer}
main{max-width:1180px;margin:0 auto;padding:12px 16px 64px}
section{margin-top:28px}
h2{font-size:19px;margin:0 0 4px;color:var(--accent);border-bottom:2px solid var(--line);padding-bottom:6px}
h2 .step{display:inline-block;background:var(--accent);color:var(--panel);border-radius:6px;padding:0 8px;margin-right:8px;font-size:14px}
h3{font-size:15px;margin:22px 0 6px}
p{max-width:900px;color:var(--ink2)} .lead{color:var(--ink2)}
.muted{color:var(--muted);font-size:12px} .bad{color:var(--bad)}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:8px;margin:10px 0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin:14px 0}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:12px 14px}
.tl{color:var(--muted);font-size:12px;text-transform:uppercase;letter-spacing:.04em}
.tv{font-size:22px;font-weight:650;margin:2px 0}
.ts{color:var(--ink2);font-size:12px}
.verdict{background:var(--panel);border:2px solid var(--accent2);border-radius:14px;padding:14px 18px;margin:10px 0}
.vt{color:var(--accent);font-weight:700;text-transform:uppercase;font-size:12px;letter-spacing:.06em}
.vm{font-size:21px;margin:2px 0} .vs{color:var(--ink2)} .vw{margin-top:6px;font-size:13px}
.warnline{margin-top:6px;color:var(--warn);font-size:13px}
.tablewrap{overflow-x:auto;margin:8px 0}
table.data{border-collapse:collapse;width:100%;font-size:13px;background:var(--panel);border:1px solid var(--line);border-radius:10px}
table.data th{position:sticky;top:0;background:var(--panel);text-align:left;color:var(--ink2);font-weight:600;border-bottom:2px solid var(--line);padding:7px 9px;white-space:nowrap}
table.data td{border-bottom:1px solid var(--line);padding:6px 9px;vertical-align:middle}
.num,table.data th.num{text-align:right;font-variant-numeric:tabular-nums}
table.grid td.rel{white-space:nowrap;font-size:12.5px}
.rel.agree{background:var(--agree-bg)} .rel.coarse{background:var(--coarse-bg)}
.rel.conflict{background:var(--conflict-bg)} .rel.nocall,.rel.nocons{background:var(--nocall-bg);color:var(--muted)}
.legend{display:flex;flex-wrap:wrap;gap:8px;font-size:12px;margin:6px 0}
.legend .key{padding:3px 8px;border-radius:6px}
.chip{display:inline-block;border-radius:999px;padding:1px 9px;font-size:12px;font-weight:600;white-space:nowrap}
.chip.good{background:var(--agree-bg);color:var(--good)} .chip.warning{background:var(--coarse-bg);color:var(--warn)}
.chip.critical{background:var(--conflict-bg);color:var(--bad)}
.meter{display:inline-block;width:60px;height:7px;border-radius:4px;background:var(--nocall-bg);vertical-align:middle;margin-right:6px;overflow:hidden}
.meter span{display:block;height:100%;background:var(--meter);border-radius:4px}
.meter.risk span{background:var(--bad)} .meter.high span{background:var(--good)}
.mv{font-variant-numeric:tabular-nums}
td:has(.meter){white-space:nowrap}
details{margin:10px 0;background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:8px 12px}
summary{cursor:pointer;font-weight:600;color:var(--ink2)}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:10px;margin-top:8px}
figure{margin:0} figure img{width:100%;border-radius:8px;background:#fff}
figcaption{font-size:12px;color:var(--muted)}
pre{font-size:12px;overflow-x:auto;background:var(--bg);padding:10px;border-radius:8px}
code{font-size:12.5px}
a{color:var(--accent)}
.refs{font-size:12.5px;color:var(--ink2);max-width:900px}
.xctl{display:flex;flex-wrap:wrap;gap:10px;align-items:center;padding:4px 6px 0}
.xctl select{font:inherit;padding:4px 6px;border-radius:8px;border:1px solid var(--line);background:var(--panel);color:var(--ink);max-width:100%}
.xctl label{font-size:13px;color:var(--ink2)}
footer{text-align:center;color:var(--muted);font-size:12px;padding:24px}
@media (max-width:700px){header h1{font-size:18px}.tv{font-size:19px}}
"""

JS = r"""
const D = window.MATCHA_EMBED;
function ink(){const s=getComputedStyle(document.documentElement);return {ink:s.getPropertyValue('--ink').trim(),ink2:s.getPropertyValue('--ink2').trim(),panel:s.getPropertyValue('--panel').trim()};}
function retheme(){const c=ink();document.querySelectorAll('.js-plotly-plot').forEach(el=>{try{Plotly.relayout(el,{'font.color':c.ink2,'legend.font.color':c.ink2});}catch(e){}});}
function buildExplorer(){
  const host=document.getElementById('explorer'); if(!host||!D) return;
  const groups={}; D.options.forEach((o,i)=>{(groups[o.group]=groups[o.group]||[]).push(i);});
  let sel='<div class="xctl"><label>Colour by <select id="xsel">';
  for(const g in groups){sel+=`<optgroup label="${g}">`+groups[g].map(i=>`<option value="${i}">${D.options[i].name}</option>`).join('')+'</optgroup>';}
  sel+='</select></label><label><input type="checkbox" id="xlab" checked> labels on the plot</label>';
  if(D.n_shown<D.n_total) sel+=`<span class="muted">${D.n_shown.toLocaleString()} of ${D.n_total.toLocaleString()} cells shown</span>`;
  sel+='</div><div id="xplot" style="height:620px"></div>';
  host.innerHTML=sel;
  const size = D.x.length>20000?3:(D.x.length>5000?4:5);
  function draw(){
    const o=D.options[+document.getElementById('xsel').value]; const showLab=document.getElementById('xlab').checked;
    const traces=[], annots=[]; const c=ink();
    if(o.kind==='cat'){
      const idx=o.cats.map(()=>[]); o.codes.forEach((k,i)=>idx[k].push(i));
      o.cats.forEach((cat,k)=>{ if(!idx[k].length) return;
        traces.push({type:'scattergl',mode:'markers',name:`${cat} (${idx[k].length.toLocaleString()})`,
          x:idx[k].map(i=>D.x[i]),y:idx[k].map(i=>D.y[i]),customdata:idx[k].map(i=>D.cluster[i]),
          marker:{size:size,color:o.colors[k],opacity:0.85},
          hovertemplate:`${cat}<br>cluster %{customdata}<extra></extra>`});});
      if(showLab && o.cats.length<=45){
        o.cats.forEach((cat,k)=>{ if(!idx[k].length||cat==='no call') return;
          const xs=idx[k].map(i=>D.x[i]).sort((a,b)=>a-b), ys=idx[k].map(i=>D.y[i]).sort((a,b)=>a-b);
          annots.push({x:xs[xs.length>>1],y:ys[ys.length>>1],text:cat,showarrow:false,
            font:{size:11,color:c.ink},bgcolor:c.panel+'cc',borderpad:2});});
      }
    } else {
      traces.push({type:'scattergl',mode:'markers',x:D.x,y:D.y,customdata:D.cluster,showlegend:false,
        marker:{size:size,color:o.values,colorscale:D.seq,colorbar:{title:{text:o.name,side:'right'},thickness:10}},
        hovertemplate:`${o.name}: %{marker.color:${o.fmt}}<br>cluster %{customdata}<extra></extra>`});
    }
    const layout={margin:{l:40,r:10,t:10,b:40},paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',
      font:{color:c.ink2,family:'system-ui, -apple-system, Segoe UI, sans-serif',size:12},
      xaxis:{title:D.basis+'1',showgrid:false,zeroline:false,showticklabels:false},
      yaxis:{title:D.basis+'2',showgrid:false,zeroline:false,showticklabels:false,scaleanchor:'x'},
      legend:{itemsizing:'constant',font:{size:11}},hovermode:'closest',annotations:annots};
    Plotly.react('xplot',traces,layout,{displaylogo:false,responsive:true});
  }
  document.getElementById('xsel').addEventListener('change',draw);
  document.getElementById('xlab').addEventListener('change',draw);
  draw();
}
function toggleTheme(){const r=document.documentElement;const dark=r.dataset.theme? r.dataset.theme==='dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  r.dataset.theme= dark?'light':'dark'; retheme(); if(document.getElementById('xsel')) document.getElementById('xsel').dispatchEvent(new Event('change'));}
window.addEventListener('load',()=>{buildExplorer();retheme();});
"""


def git_version():
    try:
        out = subprocess.run(["git", "-c", f"safe.directory={SCRIPT_DIR.parents[1]}", "-C",
                              str(SCRIPT_DIR), "describe", "--always", "--dirty", "--tags"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or "unknown"
    except Exception:  # noqa: BLE001
        return "unknown"


def build(run: Run, sample: str, settings=None):
    embed = explorer_data(run)
    sections = [
        ("overview", "Overview", None, overview(run)),
        ("step1", "QC and clustering stability", "Step 1", step1(run)),
        ("step2", "Annotators", "Step 2", step2(run)),
        ("step3", "Consensus", "Step 3", step3(run)),
        ("methods", "Methods and settings", None, methods(run, settings)),
        ("files", "Output files", None, files_section(run)),
    ]
    nav = "".join(f'<a href="#{i}">{esc(s + " · " if s else "")}{esc(t)}</a>' for i, t, s, _ in sections)
    body = "".join(
        f'<section id="{i}"><h2>{f"<span class=step>{esc(s)}</span>" if s else ""}{esc(t)}</h2>{b}</section>'
        for i, t, s, b in sections)
    when = datetime.now().strftime("%Y-%m-%d %H:%M")
    ver = git_version()
    embed_json = json.dumps(embed, separators=(",", ":")).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MatchACell report · {esc(sample)}</title>
<style>{CSS}</style>
<script>{pyo.get_plotlyjs()}</script>
</head><body>
<header><h1>MatchACell report · {esc(sample)}</h1>
<div class="sub">QC and cluster stability → annotation → consensus on the Cell Ontology ·
{len(run.obs):,} cells · generated {when} · MatchACell {esc(ver)}</div></header>
<nav>{nav}<span class="sp"></span><button onclick="toggleTheme()" title="Switch light/dark">◐ theme</button></nav>
<main>{body}</main>
<footer>Brewed with MatchACell · matcha-grade single-cell consensus</footer>
<script>window.MATCHA_EMBED={embed_json};</script>
<script>{JS}</script>
</body></html>"""


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--sample", required=True)
    p.add_argument("--step1", required=True, help="results/<sample>/matchacell")
    p.add_argument("--consensus", help="Step 3 folder (default: <step1>/consensus)")
    p.add_argument("--markers", help="Signature workbook (.xlsx) for the marker dot plot")
    p.add_argument("--settings", help="The run's configuration, as JSON or a JSON file")
    p.add_argument("-o", "--output", required=True)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    step1_dir = Path(args.step1)
    settings = None
    if args.settings:
        text = args.settings if args.settings.lstrip().startswith("{") else Path(args.settings).read_text()
        settings = json.loads(text)
    run = Run(step1_dir, Path(args.consensus) if args.consensus else step1_dir / "consensus",
              markers=args.markers or None, settings=settings)
    out = Path(args.output)
    out.write_text(build(run, args.sample, settings), encoding="utf-8")
    print(f"[report] {out} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
