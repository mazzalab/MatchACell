#!/usr/bin/env python3
"""scAdam v2 annotation with MatchACell's plots and cluster mapping outputs."""
from __future__ import annotations

import argparse
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from functions_annot import extract_best_res
from scipy import sparse

PREFIX = "scparadise_"


def prepare_query(adata, model_genes, layer):
    """Select model genes without changing the full clustered AnnData."""
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Input AnnData must contain cells and genes.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Input cell and gene names must be unique.")
    if layer and layer not in adata.layers:
        raise ValueError(f"Expression layer {layer!r} is missing; use --layer='' to use .X.")
    common = adata.var_names.intersection(model_genes, sort=False)
    if len(common) == 0:
        raise ValueError("No genes in common with the scAdam model; check species and gene identifiers.")
    fraction = len(common) / len(model_genes)
    print(f"[scParadise] {len(common)}/{len(model_genes)} model genes matched ({fraction:.1%}).")
    if fraction < 0.8:
        warnings.warn(
            "Fewer than 80% of model genes are present. Check species/gene identifiers; "
            "consider rerunning Step 1 with --use-hvg no.",
            stacklevel=2,
        )
    view = adata[:, common]
    matrix = view.layers[layer] if layer else view.X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("scAdam requires finite, non-negative log-normalized expression.")
    # scAdam densifies its input. Limit that allocation to shared model genes,
    # and leave counts, other genes, .raw and embeddings on the original object.
    return AnnData(X=matrix.astype(np.float32), obs=pd.DataFrame(index=adata.obs_names.copy()),
                   var=pd.DataFrame(index=common)), fraction


def transfer_predictions(adata, predicted, label_level, threshold):
    """Keep every annotation level and expose the selected level like scANVI."""
    if not adata.obs_names.equals(predicted.obs_names):
        raise ValueError("scAdam changed cell identities or order.")
    levels = sorted(
        int(match.group(1)) for col in predicted.obs
        if (match := re.fullmatch(r"scparadise_celltype_l(\d+)", col))
    )
    if not levels:
        raise ValueError("scAdam returned no annotation levels; use a scAdam v2 model.")
    selected = label_level or levels[-1]
    if selected not in levels:
        raise ValueError(f"Annotation level {selected} is unavailable; model returned {levels}.")
    label_col = f"{PREFIX}celltype_l{selected}"
    prob_col = f"{label_col}_probability"
    if prob_col not in predicted.obs:
        raise ValueError(f"scAdam did not return {prob_col}.")
    probabilities = pd.to_numeric(predicted.obs[prob_col], errors="raise")
    if not np.isfinite(probabilities).all() or not probabilities.between(0, 1).all():
        raise ValueError("scAdam returned invalid prediction probabilities.")
    for col in predicted.obs:
        if col.startswith(PREFIX):
            adata.obs[col] = predicted.obs[col]
        elif col in {"gradient_score", "entropy_score", "distance_score"}:
            adata.obs[f"{PREFIX}{col}"] = predicted.obs[col]
    # Honor the boolean detector mask even if upstream truncates the string
    # 'Unknown' in a numpy array whose original cell labels were shorter.
    unknown = predicted.obs.get(f"{PREFIX}unknown", pd.Series(False, index=predicted.obs_names))
    for level in levels:
        col = f"{PREFIX}celltype_l{level}"
        adata.obs[col] = adata.obs[col].astype(str).mask(unknown, "Unknown")
    adata.obs["predicted_celltype"] = adata.obs[label_col]
    adata.obs["prediction_probability"] = probabilities
    adata.obs["cell_type_pred"] = adata.obs[label_col].mask(probabilities < threshold, "Unknown")
    return selected


def cluster_mapping(adata, leiden_col):
    counts = pd.crosstab(adata.obs[leiden_col], adata.obs["cell_type_pred"])
    percentages = counts.div(counts.sum(axis=1), axis=0) * 100
    counts["dominant_cell_type"] = percentages.idxmax(axis=1)
    counts["dominant_pct"] = percentages.max(axis=1)
    return counts.rename_axis("cluster").reset_index()


def main(raw_data_file, verdict_file, model_dir, output_dir, layer="lognorm", label_level=0,
         thr=0.0, detect_unknown=False, unknown_method="voting", unknown_threshold=2.0,
         batch_size=256, device="auto", n_threads=8, seed=0):
    if not 0 <= thr <= 1 or label_level < 0:
        raise ValueError("thr must be in [0, 1] and label_level must be non-negative.")
    if batch_size < 1 or n_threads < 1 or not np.isfinite(unknown_threshold) or unknown_threshold <= 0:
        raise ValueError("batch_size, n_threads and unknown_threshold must be positive.")
    if device not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be auto, cpu or cuda.")
    if not model_dir:
        raise ValueError("Configure matchacell_annotation.scparadise.model_dir first.")
    model_path = Path(model_dir).expanduser().absolute()
    if not model_path.name.endswith("scAdam"):
        raise ValueError("The model directory name must end in scAdam (as created by download_model).")
    checkpoint_path = model_path / "model_v2.pth"
    if not checkpoint_path.is_file():
        raise FileNotFoundError(f"scAdam v2 checkpoint not found: {checkpoint_path}")
    if detect_unknown and not (model_path / "unknown_detector.json").is_file():
        raise FileNotFoundError("detect_unknown requires unknown_detector.json in the model directory.")

    # Set limits before importing the inference/plotting libraries. In 1.1.0,
    # predict initially loads with device='auto' even when CPU is requested.
    os.environ["OMP_NUM_THREADS"] = str(n_threads)
    os.environ["MKL_NUM_THREADS"] = str(n_threads)
    if device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
    import torch

    torch.set_num_threads(n_threads)
    torch.manual_seed(seed)
    np.random.seed(seed)
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA was requested but is unavailable; use device: cpu or auto.")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc
    from functions_annot import multi_umap

    adata = sc.read_h5ad(raw_data_file)
    leiden_col = extract_best_res(verdict_file)
    if leiden_col not in adata.obs or adata.obs[leiden_col].isna().any():
        raise ValueError(f"Verdict cluster column {leiden_col!r} is missing or contains missing labels.")
    # Inspect model metadata before prediction so incompatible inputs fail
    # clearly and scAdam only densifies the required genes.
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model_genes = checkpoint.get("var_names")
    if model_genes is None or len(model_genes) == 0:
        raise ValueError("The scAdam checkpoint does not contain gene names.")
    del checkpoint
    query, overlap = prepare_query(adata, model_genes, layer)

    import scparadise

    predicted = scparadise.scadam.predict(
        query, path_model=str(model_path), batch_size=batch_size, device=device,
        prefix=PREFIX, detect_unknown=detect_unknown, method=unknown_method,
        threshold=unknown_threshold,
    )
    selected = transfer_predictions(adata, predicted, label_level, thr)
    adata.uns["scparadise"] = {
        "version": scparadise.__version__, "model_dir": str(model_path),
        "layer": layer or "X", "label_level": selected, "threshold": thr,
        "detect_unknown": detect_unknown, "unknown_method": unknown_method,
        "unknown_threshold": unknown_threshold, "model_gene_overlap": overlap,
        "device": device, "batch_size": batch_size, "n_threads": n_threads, "seed": seed,
    }
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    cluster_mapping(adata, leiden_col).to_excel(
        output_path / "cluster_to_celltype_mapping.xlsx", index=False,
    )
    if "X_umap" in adata.obsm:
        sc.settings.figdir = str(output_path)
        sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)
        sc.pl.umap(adata, color="cell_type_pred", show=False,
                   save=f"_{leiden_col}_cell_type_pred.png")
        multi_umap(adata, leiden_col, "cell_type_pred", output_path)
        plt.close("all")
    else:
        print("[scParadise] No X_umap embedding; annotation and cluster table saved without plots.")
    adata.write_h5ad(output_path / "scparadise_annotated.h5ad")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="Clustered AnnData (.h5ad)")
    parser.add_argument("-v", "--verdict", required=True, help="MatchA_Verdict.txt from Step 1")
    parser.add_argument("-m", "--model_dir", required=True, help="Local scAdam v2 model folder")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--layer", default="lognorm", help="Log-normalized layer; empty uses .X")
    parser.add_argument("--label_level", type=int, default=0, help="0: finest level; 1, 2, ...: explicit level")
    parser.add_argument("-t", "--thr", type=float, default=0.0, help="Minimum label probability")
    parser.add_argument("--detect_unknown", action="store_true", help="Enable scAdam unknown detector")
    parser.add_argument("--unknown_method", default="voting",
                        choices=["voting", "gradient", "entropy", "distance", "combined"])
    parser.add_argument("--unknown_threshold", type=float, default=2.0, help="Ashman separation threshold")
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--n_threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = vars(parse_args())
    args["raw_data_file"] = args.pop("input")
    args["verdict_file"] = args.pop("verdict")
    args["output_dir"] = args.pop("output")
    main(**args)
