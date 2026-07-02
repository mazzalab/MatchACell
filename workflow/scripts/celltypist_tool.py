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
    output_dir=Path(output_dir)
    output_celltypist=Path(output_dir)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)
    
    # =====================
    # CellTypist
    # =====================

    sc.settings.figdir = output_celltypist

    # Scarica/Aggiorna la lista dei modelli
    models.download_models()

    # Visualizza i modelli disponibili per scegliere il più adatto
    print(models.models_description())

    #model_name = ["Mouse_Dendritic_Subtypes.pkl", "Human_Placenta_Decidua.pkl", "Pan_Fetal_Human.pkl"]
    #['Immune_All_Low.pkl', 'Immune_All_High.pkl', 'Adult_Human_Skin.pkl']

    for model in models_list:

        print(f"Using model: {model}")
        # Esegui l'annotazione
        predictions = celltypist.annotate(adata, 
                                        model = model, 
                                        majority_voting = True,
                                        over_clustering = leiden_col)

        # Inseriamo i risultati nel nostro oggetto AnnData
        adata.obs[f'celltypist_label_majority_voting_{model}'] = predictions.predicted_labels['majority_voting']
        
        sc.pl.umap(
            adata, 
            color=[f"celltypist_label_majority_voting_{model}"],  
            legend_fontsize=6,
            save=f"_{leiden_col}_celltypist_label_majority_voting_{model}.png")        

        multi_umap(adata, leiden_col, f'celltypist_label_majority_voting_{model}', output_dir)
        # 5. Salva il file in formato Excel (.xlsx)
        # output_unified_file = output_dir / "cluster_to_celltype_mapping.xlsx"
        # df_unified.to_excel(output_unified_file, index=False)

        # =====================
        # Save AnnData with scores
        # =====================
        output_h5ad = output_dir / "celltypist_annotated.h5ad"
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