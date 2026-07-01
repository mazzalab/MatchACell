#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

from functions_annot import rearrange_sign, extract_best_res


def main(raw_data_file, annot_file, thr, verdict_file, output_dir):

    # =====================
    # Extract best res
    # =====================
    leiden_col=extract_best_res(verdict_file)

    # =====================
    # Read input
    # =====================
    raw_data_filtered = sc.read_h5ad(raw_data_file)

    signature_dict_tot = rearrange_sign(raw_data_filtered, annot_file)

    # Remove empty signatures
    signature_dict_tot = {
        k: v for k, v in signature_dict_tot.items() if len(v) > 0
    }

    clean_sig_dict = {}

    for cell_type, genes in signature_dict_tot.items():

        existing_genes = [
            g for g in genes if g in raw_data_filtered.var_names
        ]

        if existing_genes:
            clean_sig_dict[cell_type] = existing_genes

    # =====================
    # Output folder
    # =====================
    output_rankgenes = Path(output_dir)
    output_rankgenes.mkdir(parents=True, exist_ok=True)

    output_violin = output_rankgenes / "violin"
    output_violin.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(output_rankgenes)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)

    # =====================
    # Gene signature scoring
    # =====================
    for signature_name, genes in clean_sig_dict.items():

        print(f"Scoring signature: {signature_name}")

        sc.tl.score_genes(
            raw_data_filtered,
            gene_list=genes,
            score_name=signature_name
        )

        global_median = raw_data_filtered.obs[signature_name].median()

        fig = sc.pl.violin(
            raw_data_filtered,
            keys=signature_name,
            groupby=leiden_col,
            stripplot=False,
            show=False
        )

        ax = fig.axes

        ax.axhline(
            global_median,
            color="red",
            linestyle="--",
            label=f"Global median ({global_median:.2f})"
        )

        ax.set_title(signature_name)
        ax.set_ylabel("Gene score")
        ax.set_ylim(-3, 8)

        plt.xticks(rotation=15, ha="right", fontsize=10)
        plt.tight_layout()

        plt.savefig(
            output_violin / f"{signature_name}_violin.png",
            dpi=300
        )
        plt.close()

    # =====================
    # Heatmap
    # =====================
    score_names = list(clean_sig_dict.keys())

    df_scores = (
        raw_data_filtered.obs
        .groupby(leiden_col)[score_names]
        .mean()
        .round(3)
    )

    plt.figure(figsize=(16, 10))

    sns.set_context("paper", font_scale=1.2)
    sns.set_style("white")

    sns.heatmap(
        df_scores,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        linewidths=0,
        annot_kws={"size": 8},
        cbar_kws={"label": "Signature Score"}
    )

    plt.title("Validation of Cell Type Identity via Signature Scores")
    plt.xlabel("Cell Type Signatures")
    plt.ylabel("Leiden Clusters")
    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(
        output_rankgenes / "Heatmap_Signature_Validation.png",
        dpi=300
    )
    plt.close()

    # =====================
    # Dotplot
    # =====================
    sc.pl.dotplot(
        raw_data_filtered,
        var_names=score_names,
        groupby=leiden_col,
        standard_scale="var",
        cmap="Blues",
        title="Signature Scores per Cluster",
        save="Rankgenes_dotplot.png"
    )

    # =====================
    # UMAP
    # =====================
    sc.pl.umap(
        raw_data_filtered,
        color=[leiden_col] + score_names,
        legend_loc="on data",
        cmap="magma",
        ncols=3,
        save="_signatures.png"
    )

    # =====================
    # Mapping automatico e Creazione Tabella di Report
    # =====================
    print("Creating annotation summary table with safety thresholds...")
    
    final_annotations = []
    highest_scores = []

    # Iteriamo riga per riga (cluster per cluster) sul DataFrame degli score
    for cluster_id, row in df_scores[score_names].iterrows():
        # 1. Seleziona solo le signature che superano la soglia (es. >= 0.1)
        valid_signatures = row[row >= float(thr)]
        
        # 2. Se ci sono signature valide, ordinale dal punteggio più alto al più basso
        if not valid_signatures.empty:
            valid_signatures_sorted = valid_signatures.sort_values(ascending=False)
            
            # Uniamo i nomi delle signature approvate separandoli da una virgola
            annotation_string = ", ".join(valid_signatures_sorted.index)
            final_annotations.append(annotation_string)
            
            # Registriamo lo score più alto in assoluto per il report
            highest_scores.append(valid_signatures_sorted.iloc[0])
        else:
            # Se nessuna signature supera la soglia, il cluster è Unknown
            final_annotations.append("Unknown")
            highest_scores.append(row.max()) # Mantiene comunque il valore massimo reale per diagnostica

    # Aggiungiamo le colonne generate al DataFrame dei punteggi
    df_scores["final_annotation"] = final_annotations
    df_scores["highest_score"] = highest_scores

    # 2. Ripristina l'indice del cluster come una colonna vera e propria nell'output
    df_report = df_scores.reset_index().rename(columns={leiden_col: "cluster"})

    # 3. Salva la tabella dei punteggi e delle annotazioni in Excel nella cartella radice dell'output
    df_report.to_excel(output_rankgenes / "cluster_annotation_summary.xlsx", index=False)
    
    print(f"Summary table saved to {output_rankgenes}/cluster_annotation_summary.xlsx")

    # 4. Applica formalmente la mappatura all'oggetto AnnData (.obs) per le analisi a valle e i grafici
    mapping_dict = df_scores["final_annotation"].to_dict()
    
    annotation_col = "cell_type_pred"
    raw_data_filtered.obs[annotation_col] = (
        raw_data_filtered.obs[leiden_col]
        .astype(str)
        .map(mapping_dict)
        .astype("category")
    )

    # =====================
    # Save AnnData with scores
    # =====================
    raw_data_filtered.write(
        output_rankgenes / "score_genes_annotated.h5ad"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate gene signature scores and generate plots."
    )

    parser.add_argument("-i", "--input", required=True, 
                        help="Input AnnData (.h5ad)"
    )

    parser.add_argument("-a", "--annot",required=True, 
                        help="Excel file containing gene signatures"
    )

    parser.add_argument("-t", "--thr",required=True, 
                        help="Threshold on gene score"
    )

    parser.add_argument("-v", "--verdict", required=True,
                        help="Verdict file"
    )

    parser.add_argument("-o", "--output", required=True, 
                        help="Output directory"
    )

    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        annot_file=args.annot,
        thr=args.thr,
        verdict_file=args.verdict,
        output_dir=args.output,
    )