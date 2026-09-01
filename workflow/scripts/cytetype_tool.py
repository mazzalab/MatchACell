#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cytetype import CyteType

import scanpy as sc

# Import dei moduli custom del workflow
from functions_annot import extract_best_res, multi_umap


def main(raw_data_file, verdict_file, output_dir, study_context, title, run_label, n_top_genes, auth_token):

    # =====================
    # Extract best res
    # =====================
    leiden_col = extract_best_res(verdict_file)

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
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)

    # =====================
    # CyteType
    # =====================
    print("[CyteType] Running CyteType Classification...")

    rank_key = f"rank_genes_{leiden_col}"
    # use_raw=False: adata.raw below holds the counts layer; rank on the
    # already log-normalized .X instead (CyteType requires log1p data).
    sc.tl.rank_genes_groups(adata, groupby=leiden_col, method="wilcoxon", key_added=rank_key, use_raw=False)
    adata.var["gene_symbols"] = adata.var_names

    # CLI flag > CYTETYPE_API_TOKEN env var > credentials stored by `cytetype setup`.
    auth_token = auth_token or os.environ.get("CYTETYPE_API_TOKEN") or None

    annotator = CyteType(
        adata,
        group_key=leiden_col,
        rank_key=rank_key,
        n_top_genes=n_top_genes,
        auth_token=auth_token,
    )

    query_filename = output_dir / "query.json"

    adata = annotator.run(
        study_context=study_context,
        metadata={"title": title, "run_label": run_label},
        save_query=True,
        query_filename=query_filename,
        show_progress=True,
    )

    # =====================
    # Mappatura dei Cluster
    # =====================
    # annotator.run() already wrote per-cell columns onto adata.obs, correctly
    # aligned to leiden_col (no need to re-derive cluster ids from query.json).
    print("CyteType Generating unified cluster annotation and counts file...")

    label_col = f"cytetype_annotation_{leiden_col}"
    report_cols = [
        c for c in (
            label_col,
            f"cytetype_cellOntologyTerm_{leiden_col}",
            f"cytetype_cellOntologyTermID_{leiden_col}",
            f"cytetype_cellState_{leiden_col}",
        )
        if c in adata.obs.columns
    ]

    df_report = (
        adata.obs.groupby(leiden_col)[report_cols]
        .first()
        .reset_index()
        .rename(columns={leiden_col: "cluster"})
    )

    output_unified_file = output_dir / "cluster_to_celltype_mapping.xlsx"
    df_report.to_excel(output_unified_file, index=False)

    sc.pl.umap(
        adata,
        color=label_col,
        save=f"_{leiden_col}_cytetype_results.png",
    )

    multi_umap(adata, leiden_col, label_col, output_dir)

    # =====================
    # Save AnnData with scores
    # =====================
    output_h5ad = output_dir / "cytetype_annotated.h5ad"
    adata.write(output_h5ad)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Launch CyteType analysis and generate plots."
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input AnnData (.h5ad)"
    )

    parser.add_argument("-v", "--verdict", required=True,
                        help="Verdict file"
    )

    parser.add_argument("-o", "--output", required=True,
                        help="Output directory"
    )

    parser.add_argument("-c", "--study_context", required=True,
                        help="Free-text study context passed to CyteType"
    )

    parser.add_argument("--title", default="MatchACell run",
                        help="Run title recorded in CyteType metadata"
    )

    parser.add_argument("--run_label", default="v1",
                        help="Run label recorded in CyteType metadata"
    )

    parser.add_argument("-n", "--n_top_genes", type=int, default=100,
                        help="Number of top marker genes passed to CyteType"
    )

    parser.add_argument("--auth_token", default=None,
                        help="CyteType API bearer token. If omitted, falls back to "
                             "credentials stored by `cytetype setup`."
    )

    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        verdict_file=args.verdict,
        output_dir=args.output,
        study_context=args.study_context,
        title=args.title,
        run_label=args.run_label,
        n_top_genes=args.n_top_genes,
        auth_token=args.auth_token,
    )
