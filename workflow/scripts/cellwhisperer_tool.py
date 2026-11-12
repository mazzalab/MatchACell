#!/usr/bin/env python3
"""
MatchACell - Step 2: CellWhisperer annotation.

CellWhisperer (Schaefer et al., Nat Biotechnol 2025) embeds transcriptomes and
free text into one space. Each cell is scored against every candidate label
(the text "{label}" filled into a template, e.g. "{label} from human pancreas")
and gets the label with the highest score. Both embeddings run on the local
CellWhisperer model, so no data leaves the machine.

This script runs inside the CellWhisperer pixi environment (see
docs/annotators/cellwhisperer.md), not in a MatchACell conda environment.
"""

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F

from functions_annot import extract_best_res, multi_umap

try:  # Excel I/O; not in CellWhisperer's own environment by default
    import openpyxl  # noqa: F401
except ImportError as e:
    raise SystemExit(
        "openpyxl is missing from the CellWhisperer environment. "
        "Run once in the CellWhisperer checkout: pixi add openpyxl"
    ) from e

# Raw counts are needed: CellWhisperer's Geneformer processor ranks genes by
# normalized expression itself.
COUNTS_LAYER = "counts"
# Cells tokenized per chunk; bounds memory on large samples.
CHUNK_CELLS = 2000


def candidate_labels(labels, annot_file):
    """Candidate labels from --labels, else the sheet names of the workbook."""
    if labels:
        out = [str(x).strip() for x in labels if str(x).strip()]
    elif annot_file:
        out = [s.strip() for s in pd.ExcelFile(annot_file).sheet_names if s.strip()]
    else:
        out = []
    out = list(dict.fromkeys(out))  # keep order, drop duplicates
    if len(out) < 2:
        raise ValueError(
            "CellWhisperer needs at least 2 candidate labels: set "
            "matchacell_annotation.cellwhisperer.labels or annot_file."
        )
    return out


def natural_key(x):
    return [(0, int(p)) if p.isdigit() else (1, p) for p in str(x).split("_")]


def embed_cells(adata, model, processor, batch_size):
    from cellwhisperer.utils.processing import adata_to_embeds

    chunks = []
    for start in range(0, adata.n_obs, CHUNK_CELLS):
        sub = adata[start:start + CHUNK_CELLS].copy()
        with torch.no_grad():
            chunks.append(adata_to_embeds(sub, model, processor, batch_size=batch_size).cpu())
        print(f"  embedded {min(start + CHUNK_CELLS, adata.n_obs)}/{adata.n_obs} cells")
    return torch.cat(chunks)


def main(args):
    torch.set_num_threads(args.n_threads)
    leiden_col = extract_best_res(args.verdict)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out

    labels = candidate_labels(args.labels, args.annot_file)
    texts = [args.template.format(label=lab) for lab in labels]
    print(f"CellWhisperer: {len(labels)} candidate labels, template {args.template!r}")

    adata = sc.read_h5ad(args.input)
    if leiden_col not in adata.obs:
        raise KeyError(f"Verdict column {leiden_col!r} not in {args.input}")
    if COUNTS_LAYER not in adata.layers:
        raise KeyError(f"Layer {COUNTS_LAYER!r} with raw counts not in {args.input}")
    if adata.n_vars < 5000:
        print(f"WARNING: only {adata.n_vars} genes. CellWhisperer ranks genes per cell, "
              "so a reduced gene set (HVGs, a targeted panel) weakens its embeddings. "
              "For Step 1 output, set matchacell_cluster_stability.extra: '--use-hvg no'.")

    query = adata.copy()
    query.X = query.layers[COUNTS_LAYER]
    if "gene_name" not in query.var:
        query.var["gene_name"] = query.var_names.astype(str)

    # cellwhisperer.config reads ./config.yaml at import time, which would pick
    # up MatchACell's config: import it from inside the CellWhisperer checkout.
    os.chdir(args.repo_dir)
    from cellwhisperer.utils.model_io import load_cellwhisperer_model

    pl_model, _, processor = load_cellwhisperer_model(args.model)
    model = pl_model.model
    if args.device != "auto":
        model.to(args.device)
    logit_scale = float(model.discriminator.temperature.exp())
    print(f"Model on {model.device}, logit scale {logit_scale:.2f}")

    with torch.no_grad():
        text_embeds = model.embed_texts(texts, chunk_size=128).float().cpu()
    cell_embeds = embed_cells(query, model, processor, args.batch_size).float()

    # Per cell: softmax over the candidates of the CLIP logits.
    logits = logit_scale * cell_embeds @ text_embeds.T          # cells x labels
    prob = logits.softmax(dim=1)
    best_p, best_i = prob.max(dim=1)
    pred = np.asarray(labels, dtype=object)[best_i.numpy()]
    pred[best_p.numpy() < args.threshold] = "Unknown"

    adata.obs["cell_type_pred"] = pd.Categorical(pred)
    adata.obs["prediction_probability"] = best_p.numpy()
    adata.obsm["X_cellwhisperer"] = cell_embeds.numpy()
    adata.uns["cellwhisperer"] = {
        "labels": labels, "template": args.template, "model": str(args.model),
        "leiden_col": leiden_col, "threshold": args.threshold,
    }
    score_df = pd.DataFrame(logits.numpy(), index=adata.obs_names,
                            columns=[f"cellwhisperer_score_{lab}" for lab in labels])

    # Per cluster: majority of the cell calls, and the call of the cluster's
    # mean embedding (one profile per cluster, like the cluster-level tools).
    rows = []
    for cl in sorted(adata.obs[leiden_col].unique(), key=natural_key):
        mask = (adata.obs[leiden_col] == cl).to_numpy()
        counts = pd.Series(pred[mask]).value_counts()
        centroid = F.normalize(cell_embeds[mask].mean(0, keepdim=True), dim=1)
        c_prob = (logit_scale * centroid @ text_embeds.T).softmax(dim=1)[0]
        rows.append({
            "cluster": cl,
            "n_cells": int(mask.sum()),
            "majority_label": counts.index[0],
            "majority_fraction": round(counts.iloc[0] / mask.sum(), 3),
            "centroid_label": labels[int(c_prob.argmax())],
            "centroid_probability": round(float(c_prob.max()), 3),
            "mean_probability": round(float(best_p[mask].mean()), 3),
        })
    table = pd.DataFrame(rows)
    table.to_excel(out / "cluster_to_celltype_mapping.xlsx", index=False)
    print(table.to_string(index=False))

    if any(k.lower() == "x_umap" for k in adata.obsm):
        sc.pl.umap(adata, color="cell_type_pred", legend_fontsize=6, show=False,
                   save=f"_{leiden_col}_cellwhisperer.png")
        multi_umap(adata, leiden_col, "cell_type_pred", out)

    adata.write_h5ad(out / "cellwhisperer_annotated.h5ad")
    score_df.to_csv(out / "cellwhisperer_cell_scores.csv.gz")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="CellWhisperer annotation for MatchACell.")
    p.add_argument("-i", "--input", required=True, help="Step 1 clustered .h5ad")
    p.add_argument("-v", "--verdict", required=True, help="MatchA_Verdict.txt")
    p.add_argument("-o", "--output", required=True, help="Output directory")
    p.add_argument("-m", "--model", required=True, help="cellwhisperer_clip_v1.ckpt")
    p.add_argument("-r", "--repo_dir", required=True, help="CellWhisperer checkout")
    p.add_argument("--labels", nargs="*", default=[], help="Candidate labels")
    p.add_argument("--annot_file", default="", help="Signature workbook; its sheet names "
                   "are the candidates when --labels is empty")
    p.add_argument("--template", default="{label}",
                   help="Text scored per label, e.g. '{label} from human pancreas'")
    p.add_argument("--threshold", type=float, default=0.0,
                   help="Minimum softmax probability; below -> Unknown")
    p.add_argument("--batch_size", type=int, default=32)
    p.add_argument("--device", default="auto", help="auto | cpu | cuda")
    p.add_argument("--n_threads", type=int, default=8)
    args = p.parse_args()
    if "{label}" not in args.template:
        p.error("--template must contain {label}")
    # Absolute paths: main() changes directory before loading the model.
    for name in ("input", "verdict", "output", "model", "repo_dir", "annot_file"):
        if getattr(args, name):
            setattr(args, name, os.path.abspath(getattr(args, name)))
    main(args)
