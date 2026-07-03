#!/usr/bin/env python3
import argparse
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cia import investigate, report, utils

import pandas as pd
import scanpy as sc

# Import dei moduli custom del workflow
from functions_annot import rearrange_sign, extract_best_res, multi_umap
# Moduli ipotizzati dallo script originale (investigate, utils, report)
# Assicurati che siano importabili nel path o aggiungi i relativi import corretti:
# import investigate
# import utils
# import report

def main(raw_data_file, annot_file, verdict_file, thr, output_dir, ncpus):

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
        
    signature_dict_tot = rearrange_sign(adata, annot_file)

    # Remove empty signatures
    signature_dict_tot = {
        k: v for k, v in signature_dict_tot.items() if len(v) > 0
    }

    clean_sig_dict = {}

    for cell_type, genes in signature_dict_tot.items():

        existing_genes = [
            g for g in genes if g in adata.var_names
        ]

        if existing_genes:
            clean_sig_dict[cell_type] = existing_genes
    
    # =====================
    # Output folder
    # =====================
    output_dir=Path(output_dir)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)
    
    # =====================
    # CIA
    # =====================
    gmt_orig = investigate.load_signatures(signature_dict_tot)

    for sig_name, genes in gmt_orig.items():
        print(f"  - {sig_name}: {len(genes)} genes")

    # Controllo di similarità delle signature (Jaccard Index)
    utils.signatures_similarity(gmt_orig, show='J')

    # Filtraggio delle signature con meno di 3 geni
    gmt = {k: v for k, v in gmt_orig.items() if len(v) >= 3}

    # 5. Classificazione CIA
    print("[CIA] Running CIA Classification...")
    investigate.CIA_classify(
        data=adata, 
        signatures_input=gmt, 
        label_column='CIA prediction default', 
        n_cpus=ncpus
    )

    # Salvataggio dei punteggi delle cellule
    df_scores = adata.obs["CIA prediction default"]
    df_scores.to_csv(output_dir / "signature_scores_cell.csv")

    # Conversione a categoria della classificazione
    adata.obs['CIA prediction default'] = adata.obs['CIA prediction default'].astype('category')
    categories = adata.obs['CIA prediction default'].cat.categories.tolist()

    # 7. Configurazione Palette Colori Dedicata
    base_palette = sc.pl.palettes.default_102
    n_colors_needed = len([c for c in categories if c != 'Unassigned'])
    assert n_colors_needed <= len(base_palette), "Not enough unique colors in base_palette"

    cia_colors = {}
    color_index = 0
    for cat in categories:
        if cat == 'Unassigned':
            cia_colors[cat] = '#cccccc'  # Grigio per le cellule non assegnate
        else:
            cia_colors[cat] = base_palette[color_index]
            color_index += 1

    color_list = [cia_colors[cat] for cat in categories]

    # 8. Generazione dei Plot (Tutti con show=False e salvati a 150 DPI)
    print("[CIA] Plotting UMAPs and metrics...")
    
    # UMAP standard con palette dedicata
    sc.pl.umap(
        adata,
        color=['CIA prediction default'],
        palette=color_list,
        save="_cia_clusters_cell.png",
        show=False
    )

    # Multi UMAP (se la funzione custom è definita nel tuo ambiente)
    multi_umap(adata, leiden_col, 'CIA prediction default', output_dir)

    # Metrics di classificazione e Group Composition
    # report.group_composition(
    #     adata, 
    #     classification_obs='CIA prediction default', 
    #     ref_obs=leiden_col,
    #     cmap='Greens',
    #     save="_group_composition.png"
    # )

    # =====================
    # Mappatura dei Cluster
    # =====================
    print("CIA Generating unified cluster annotation and counts file...")
    
    # 1. Crea la tabella di contingenza standard (conteggi assoluti)
    df_unified = pd.crosstab(adata.obs[leiden_col], adata.obs['CIA prediction default'])
    
    # 2. Calcola le percentuali per applicare la soglia dell'80%
    percentages = df_unified.div(df_unified.sum(axis=1), axis=0) * 100

    cluster_annotations = []
    for cluster in percentages.index:
        row_pct = percentages.loc[cluster]
        
        # Trova tutte le signature che superano o uguagliano la soglia del 70%
        # Escludiamo 'Unassigned' dal potenziale nome del tipo cellulare
        valid_sigs = [
            sig for sig, pct in row_pct.items() 
            if pct >= float(thr) and sig != 'Unassigned'
        ]

        # Se ci sono signature valide le uniamo con una virgola, altrimenti "Ambiguous/Mixed"
        if len(valid_sigs) > 0:
            cluster_annotations.append(", ".join(valid_sigs))
        else:
            cluster_annotations.append("Unknown")

    # 3. Aggiunge la colonna con il verdetto ad alta flessibilità
    df_unified['CIA_cluster_annotation'] = cluster_annotations

    # 4. Sposta il cluster dall'indice a una colonna vera e propria
    df_unified = df_unified.reset_index()
    df_unified = df_unified.rename(columns={leiden_col: 'cluster'})

    # 5. Salva il file in formato Excel (.xlsx)
    output_unified_file = output_dir / "cluster_to_celltype_mapping.xlsx"
    df_unified.to_excel(output_unified_file, index=False)

    # =====================
    # Save AnnData with scores
    # =====================
    output_h5ad = output_dir / "CIA_annotated.h5ad"
    adata.write(output_h5ad)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Launch CIA analysis and generate plots."
    )

    parser.add_argument("-i", "--input", required=True, 
                        help="Input AnnData (.h5ad)"
    )

    parser.add_argument("-a", "--annot",required=True, 
                        help="Excel file containing gene signatures"
    )

    parser.add_argument("-v", "--verdict", required=True,
                        help="Verdict file"
    )

    parser.add_argument("-o", "--output", required=True, 
                        help="Output directory"
    )

    parser.add_argument("-t", "--thr", required=True, 
                        help="Threshold on CIA"
    )

    parser.add_argument("-n", "--ncpus", required=True, type=int, 
                        help="Number of CPUs"
    )

    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        annot_file=args.annot,
        verdict_file=args.verdict,
        output_dir=args.output,
        thr=args.thr,
        ncpus=args.ncpus,
    )