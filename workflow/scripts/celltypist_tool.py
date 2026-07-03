#!/usr/bin/env python3
import argparse
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import celltypist
from celltypist import models

import pandas as pd
import scanpy as sc

# Import dei moduli custom del workflow
from functions_annot import extract_best_res, multi_umap
# Moduli ipotizzati dallo script originale (investigate, utils, report)
# Assicurati che siano importabili nel path o aggiungi i relativi import corretti:
# import investigate
# import utils
# import report

def main(raw_data_file, verdict_file, models_list, output_dir):

    # =====================
    # Extract best res
    # =====================
    leiden_col=extract_best_res(verdict_file)

    # =====================
    # Read input
    # =====================
    adata = sc.read_h5ad(raw_data_file)
    adata_counts = adata.copy()
    adata_counts.X = adata.layers["counts"]
    adata.raw = adata_counts
    
    # =====================
    # Output folder
    # =====================
    output_celltypist=Path(output_dir)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)
    
    # =====================
    # CellTypist
    # =====================

    sc.settings.figdir = output_celltypist

    # Scarica/Aggiorna la lista dei modelli
    models.download_models()

    unique_clusters = sorted(adata.obs[leiden_col].unique(), key=lambda x: int(x) if str(x).isdigit() else x)
    df_final_report = pd.DataFrame({'cluster': unique_clusters})
    df_final_report = df_final_report.set_index('cluster')

    for model in models_list:

        print(f"Using model: {model}")
        # Esegui l'annotazione
        predictions = celltypist.annotate(adata, 
                                        model = model, 
                                        majority_voting = True,
                                        over_clustering = leiden_col)

        # Inseriamo i risultati nel nostro oggetto AnnData
        label_col_name = f'celltypist_label_majority_voting_{model}'
        adata.obs[label_col_name] = predictions.predicted_labels['majority_voting']

        sc.pl.umap(
            adata, 
            color=label_col_name,  
            legend_fontsize=6,
            save=f"_{leiden_col}_{label_col_name}.png")        

        multi_umap(adata, leiden_col, label_col_name, output_celltypist)

        col_name_excel = model.replace(".pkl", "")
        majority_labels = predictions.predicted_labels['majority_voting']
        cluster_mapping = adata.obs.groupby(leiden_col)[f'celltypist_label_majority_voting_{model}'].first()
        
        # Mappiamo i dati sul report finale
        col_name_excel = model.replace(".pkl", "")
        df_final_report[col_name_excel] = df_final_report.index.map(cluster_mapping)

    
    # =====================
    # 6. Salvataggio Report Excel Unico Foglio
    # =====================
    df_final_report = df_final_report.reset_index()
    output_excel_file = output_celltypist / "cluster_to_celltype_mapping.xlsx"
    df_final_report.to_excel(output_excel_file, index=False)

    # =====================
    # Save AnnData with scores
    # =====================
    output_h5ad = output_celltypist / "celltypist_annotated.h5ad"
    adata.write(output_h5ad)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Launch Celltypist analysis and generate plots."
    )

    parser.add_argument("-i", "--input", required=True, 
                        help="Input AnnData (.h5ad)"
    )

    parser.add_argument("-v", "--verdict", required=True,
                        help="Verdict file"
    )

    parser.add_argument("-m", "--models_list", required=True, nargs='+', 
                        help="celltypist models list"
    )

    parser.add_argument("-o", "--output", required=True,
                        help="Output directory"
    )
    
    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        verdict_file=args.verdict,
        models_list=args.models_list,
        output_dir=args.output,
    )
