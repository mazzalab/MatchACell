#!/usr/bin/env python3
import argparse
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import scanpy as sc
import scvi
import torch

# Import dei moduli custom del workflow
from functions_annot import extract_best_res, multi_umap


def main(raw_data_file, verdict_file, reference_file, output_dir, label_col, unlabeled_category, thr,
         max_epochs_ref, max_epochs_scanvi, max_epochs_query, batch_size, n_threads):

    # =====================
    # Extract best res
    # =====================
    leiden_col = extract_best_res(verdict_file)

    # =====================
    # Read input
    # =====================
    adata = sc.read_h5ad(raw_data_file)

    # scVI/SCANVI model raw counts (negative-binomial likelihood), not the
    # log-normalized matrix -- both are already on disk as layers (see
    # workflow/scripts/matchacell_cluster_stability.py), no need to re-derive
    # them from an external zarr store like the original script did.
    adata.layers["lognorm"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()

    # =====================
    # Output folder
    # =====================
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)

    # =====================
    # Reference atlas
    # =====================
    print(f"[scANVI] Loading reference atlas: {reference_file}")
    adata_ref = sc.read_h5ad(reference_file)

    # Reference atlases commonly keep raw counts in .raw (with normalized
    # data in .X); recover them for scVI's raw-count model, same as query.
    if adata_ref.raw is not None:
        adata_ref = adata_ref.raw.to_adata()

    # Atlases often leave unannotated cells either as a real NaN or -- as
    # Tabula Muris Senis does -- as the literal string "nan" written into
    # the category itself (isna() misses that: it's a valid category, not
    # a missing value). Left alone, "nan" trains as its own class and, if
    # it happens to be the majority label like here (55% of the reference),
    # SCANVI collapses to predicting it for nearly every query cell. Catch
    # both forms.
    raw_labels = adata_ref.obs[label_col].astype(object)
    is_missing = raw_labels.isna() | raw_labels.astype(str).str.strip().str.lower().eq("nan")
    n_missing = int(is_missing.sum())
    if n_missing:
        print(f"[scANVI] {n_missing} reference cells have no real '{label_col}' label -> '{unlabeled_category}'")
        raw_labels[is_missing] = unlabeled_category
        adata_ref.obs[label_col] = raw_labels.astype("category")

    # Keep only genes shared between reference and query
    common_genes = adata_ref.var_names.intersection(adata.var_names)
    if len(common_genes) == 0:
        raise ValueError(
            f"No genes in common between reference ({reference_file}) and query ({raw_data_file})."
        )
    print(f"[scANVI] {len(common_genes)} genes in common with the reference.")

    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata[:, common_genes].copy()

    # =====================
    # scVI -> SCANVI reference mapping
    # =====================
    torch.set_num_threads(n_threads)
    scvi.settings.num_threads = n_threads

    scvi.model.SCVI.setup_anndata(adata_ref, batch_key=None, labels_key=label_col)

    print("[scANVI] Training SCVI on the reference...")
    vae_ref = scvi.model.SCVI(adata_ref, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae_ref.train(max_epochs=max_epochs_ref, early_stopping=True, batch_size=batch_size)

    print("[scANVI] Converting to SCANVI classifier...")
    lvae_ref = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        labels_key=label_col,
        unlabeled_category=unlabeled_category,
    )
    lvae_ref.train(max_epochs=max_epochs_scanvi, batch_size=batch_size)

    print("[scANVI] Mapping query onto the reference...")
    vae_q = scvi.model.SCANVI.load_query_data(adata_query, lvae_ref)
    vae_q.train(max_epochs=max_epochs_query, plan_kwargs={"weight_decay": 0.0}, batch_size=batch_size)

    adata_query.obs["predicted_celltype"] = vae_q.predict()
    adata_query.obs["prediction_probability"] = vae_q.predict(soft=True).max(axis=1)

    # =====================
    # Trasferimento delle predizioni sull'oggetto query completo
    # =====================
    adata.obs["predicted_celltype"] = adata_query.obs["predicted_celltype"].astype(str)
    adata.obs["prediction_probability"] = adata_query.obs["prediction_probability"]

    label_col_out = "cell_type_pred"
    adata.obs[label_col_out] = adata.obs["predicted_celltype"]
    low_confidence = adata.obs["prediction_probability"] < float(thr)
    adata.obs.loc[low_confidence, label_col_out] = "Unknown"

    n_unknown = int(low_confidence.sum())
    pct_unknown = 100 * n_unknown / adata.n_obs
    print(f"[scANVI] {n_unknown} cells ({pct_unknown:.1f}%) below confidence threshold {thr} -> Unknown")

    # Ripristina i dati log-normalizzati per i plot a valle
    adata.X = adata.layers["lognorm"]

    # =====================
    # UMAP
    # =====================
    sc.pl.umap(
        adata,
        color=label_col_out,
        save=f"_{leiden_col}_{label_col_out}.png",
    )

    multi_umap(adata, leiden_col, label_col_out, output_dir)

    # =====================
    # Mappatura dei Cluster
    # =====================
    print("scANVI Generating unified cluster annotation and counts file...")

    df_unified = pd.crosstab(adata.obs[leiden_col], adata.obs[label_col_out])
    percentages = df_unified.div(df_unified.sum(axis=1), axis=0) * 100

    dominant_celltype = []
    dominant_pct = []
    for cluster in percentages.index:
        row_pct = percentages.loc[cluster]
        top_label = row_pct.idxmax()
        dominant_celltype.append(top_label)
        dominant_pct.append(row_pct[top_label])

    df_unified["dominant_cell_type"] = dominant_celltype
    df_unified["dominant_pct"] = dominant_pct
    df_unified = df_unified.reset_index().rename(columns={leiden_col: "cluster"})

    output_unified_file = output_dir / "cluster_to_celltype_mapping.xlsx"
    df_unified.to_excel(output_unified_file, index=False)

    # =====================
    # Save AnnData with predictions
    # =====================
    output_h5ad = output_dir / "scanvi_annotated.h5ad"
    adata.write(output_h5ad)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Launch scANVI reference-mapping annotation and generate plots."
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input AnnData (.h5ad)"
    )

    parser.add_argument("-v", "--verdict", required=True,
                        help="Verdict file"
    )

    parser.add_argument("-r", "--reference", required=True,
                        help="Reference annotated AnnData (.h5ad) used to train the scVI/SCANVI model"
    )

    parser.add_argument("-o", "--output", required=True,
                        help="Output directory"
    )

    parser.add_argument("--label_col", default="free_annotation",
                        help="Column in the reference .obs holding cell type labels"
    )

    parser.add_argument("--unlabeled_category", default="Unknown",
                        help="Placeholder label SCANVI uses for unlabeled cells"
    )

    parser.add_argument("-t", "--thr", type=float, default=0.5,
                        help="Minimum prediction probability to keep a label (below -> Unknown)"
    )

    parser.add_argument("--max_epochs_ref", type=int, default=20,
                        help="Max training epochs for the reference SCVI model"
    )

    parser.add_argument("--max_epochs_scanvi", type=int, default=20,
                        help="Max training epochs for the SCANVI classifier"
    )

    parser.add_argument("--max_epochs_query", type=int, default=10,
                        help="Max training epochs for query-to-reference mapping"
    )

    parser.add_argument("--batch_size", type=int, default=2048,
                        help="Training batch size"
    )

    parser.add_argument("--n_threads", type=int, default=10,
                        help="Number of CPU threads for torch/scvi"
    )

    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        verdict_file=args.verdict,
        reference_file=args.reference,
        output_dir=args.output,
        label_col=args.label_col,
        unlabeled_category=args.unlabeled_category,
        thr=args.thr,
        max_epochs_ref=args.max_epochs_ref,
        max_epochs_scanvi=args.max_epochs_scanvi,
        max_epochs_query=args.max_epochs_query,
        batch_size=args.batch_size,
        n_threads=args.n_threads,
    )
