#!/usr/bin/env python3
import argparse
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import anndata as ad
import numpy as np
import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import panhumanpy as ph

# Import dei moduli custom del workflow (functions_annot deve essere nel PYTHONPATH);
# se manca, lo script gira comunque in modalità standalone con --cluster_col
try:
    from functions_annot import extract_best_res, multi_umap
except ModuleNotFoundError:
    extract_best_res = multi_umap = None

PRED_COL = 'Azimuth prediction'
UNASSIGNED_TOKENS = {'False', 'nan', 'None', ''}


def main(raw_data_file, verdict_file, cluster_col, level, conf_thr, thr, output_dir, model_version):

    # =====================
    # Extract best res
    # =====================
    if cluster_col:
        leiden_col = cluster_col
    elif verdict_file and extract_best_res is not None:
        leiden_col = extract_best_res(verdict_file)
    else:
        raise SystemExit("Serve --verdict (con functions_annot importabile) oppure --cluster_col")

    # =====================
    # Read input
    # =====================
    adata = sc.read_h5ad(raw_data_file)

    # Azimuth gira su un oggetto minimale con i raw counts (tutti i geni);
    # adata originale (X normalizzata, UMAP, leiden) resta intatto per plot e mapping
    if "counts" in adata.layers:
        counts = adata.layers["counts"]
    else:
        print("[Azimuth] layers['counts'] non trovato: uso .X come counts")
        counts = adata.X
    vals = counts.data[:100_000] if sp.issparse(counts) else np.asarray(counts).ravel()[:100_000]
    if (vals < 0).any() or not np.allclose(vals, np.round(vals)):
        raise SystemExit("I counts usati non sono interi non negativi: servono raw counts "
                         "in layers['counts'] (o in .X)")

    adata_counts = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    adata_counts.var_names_make_unique()

    # =====================
    # Output folder
    # =====================
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)

    # =====================
    # Azimuth
    # =====================
    print("[Azimuth] Running Azimuth annotation...")
    kw = {"model_version": model_version} if model_version else {}
    azimuth = ph.AzimuthNN(adata_counts, **kw)

    for col in ("azimuth_broad", "azimuth_medium", "azimuth_fine"):
        azimuth.map_to_cell_ontology(col, include_cl_id=True)

    meta = azimuth.cells_meta
    new_cols = [c for c in meta.columns if c not in adata.obs.columns]
    adata.obs = adata.obs.join(meta[new_cols])

    # Predizione finale: livello scelto, 'False'/NaN e bassa confidenza -> Unassigned
    pred = adata.obs[level].astype(str)
    pred[pred.isin(UNASSIGNED_TOKENS)] = 'Unassigned'
    if conf_thr is not None:
        low = adata.obs['final_level_confidence'].astype(float) < conf_thr
        pred[low] = 'Unassigned'
        print(f"[Azimuth] {low.sum()} cells below confidence {conf_thr} -> Unassigned")
    adata.obs[PRED_COL] = pred.astype('category')

    # Salvataggio delle predizioni per cellula (tutte le colonne Azimuth + predizione finale)
    adata.obs[new_cols + [PRED_COL]].to_csv(output_dir / "azimuth_prediction_cell.csv")

  
    categories = adata.obs[PRED_COL].cat.categories.tolist()
    base_palette = sc.pl.palettes.default_102
    n_colors_needed = len([c for c in categories if c != 'Unassigned'])
    assert n_colors_needed <= len(base_palette), "Not enough unique colors in base_palette"

    colors, i = [], 0
    for cat in categories:
        if cat == 'Unassigned':
            colors.append('#cccccc')
        else:
            colors.append(base_palette[i])
            i += 1

    # =====================
    # Plot
    # =====================
    print("[Azimuth] Plotting UMAPs...")
    sc.pl.umap(
        adata,
        color=[PRED_COL],
        palette=colors,
        save="_azimuth_clusters_cell.png",
        show=False
    )
    if multi_umap is not None:
        multi_umap(adata, leiden_col, PRED_COL, output_dir)
    else:
        sc.pl.umap(adata, color=[leiden_col, PRED_COL], wspace=0.4,
                   save="_azimuth_vs_clusters.png", show=False)

    # =====================
    # Cluster Mapping
    # =====================
    print("[Azimuth] Generating unified cluster annotation and counts file...")
    df_unified = pd.crosstab(adata.obs[leiden_col], adata.obs[PRED_COL])
    percentages = df_unified.div(df_unified.sum(axis=1), axis=0) * 100

    cluster_annotations = []
    for cluster in percentages.index:
        row_pct = percentages.loc[cluster]
        valid = [ct for ct, pct in row_pct.items() if pct >= thr and ct != 'Unassigned']
        cluster_annotations.append(", ".join(valid) if valid else "Unknown")

    df_unified['Azimuth_cluster_annotation'] = cluster_annotations
    df_unified = df_unified.reset_index().rename(columns={leiden_col: 'cluster'})
    df_unified.to_excel(output_dir / "cluster_to_celltype_mapping.xlsx", index=False)

    # =====================
    # Save AnnData
    # =====================
    adata.write(output_dir / "Azimuth_annotated.h5ad")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Launch Azimuth (panhumanpy) annotation and generate plots."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input AnnData (.h5ad) with raw counts in layers['counts']")
    parser.add_argument("-v", "--verdict", default=None,
                        help="Verdict file (richiede functions_annot)")
    parser.add_argument("--cluster_col", default=None,
                        help="Colonna .obs dei cluster da usare al posto del verdict (es. louvain)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")
    parser.add_argument("-l", "--level", default="azimuth_fine",
                        choices=["azimuth_broad", "azimuth_medium", "azimuth_fine", "final_level_labels"],
                        help="Azimuth annotation level used as prediction")
    parser.add_argument("-c", "--conf", type=float, default=None,
                        help="Optional min final_level_confidence (cells below -> Unassigned). "
                             "Default: no filter, confidence is exported per cell")
    parser.add_argument("-t", "--thr", required=True, type=float,
                        help="Min %% of cluster cells for cluster annotation")
    parser.add_argument("-m", "--model_version", default=None,
                        help="Azimuth model version (default: package default)")

    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        verdict_file=args.verdict,
        cluster_col=args.cluster_col,
        level=args.level,
        conf_thr=args.conf,
        thr=args.thr,
        output_dir=args.output,
        model_version=args.model_version,
    )