#!/usr/bin/env python3
import argparse
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import celltypeai as cta

import scanpy as sc

# Import dei moduli custom del workflow
from functions_annot import extract_best_res, multi_umap


def main(raw_data_file, verdict_file, output_dir, species, tissue, model, num_genes, n_iterations, verbose):

    # =====================
    # Extract best res
    # =====================
    leiden_col = extract_best_res(verdict_file)

    # =====================
    # Read input
    # =====================
    adata = sc.read_h5ad(raw_data_file)

    # Unlike celltypist/cytetype, this tool never touches adata.raw: its
    # internal sc.tl.rank_genes_groups() call (in cta.cell_annotator) does
    # not expose use_raw=False, and sc.tl.rank_genes_groups() defaults to
    # adata.raw when it is set -- ranking on raw counts instead of the
    # log-normalized .X already on disk (and triggering scanpy's "rank_genes_
    # groups on raw count data" warning). Leaving adata.raw unset makes it
    # rank on .X, which is what we want.

    # =====================
    # Output folder
    # =====================
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150, dpi_save=150, frameon=False)

    # =====================
    # CellTypeAI
    # =====================
    print("[CellTypeAI] Running CellTypeAI ensemble annotation via local Ollama...")

    # cell_annotator() only ever looks for a column literally named "leiden"
    # or "louvain" (github.com/rhdaw/CellTypeAI), so alias the workflow's
    # chosen resolution onto "leiden" regardless of what MatchA_Verdict.txt
    # actually recommended. It also runs its own internal
    # sc.tl.rank_genes_groups(groupby="leiden", method="wilcoxon") on
    # adata.X, so this needs to already be log-normalized (as elsewhere in
    # the workflow) rather than raw counts.
    adata.obs["leiden"] = adata.obs[leiden_col].astype(str)

    adata = cta.cell_annotator(
        species,
        tissue,
        adata,
        model=model,
        num_genes=num_genes,
        n_iterations=n_iterations,
        verbose=verbose,
    )

    # =====================
    # Mappatura dei Cluster
    # =====================
    print("CellTypeAI Generating unified cluster annotation and counts file...")

    label_col = "cell_type_ai"

    df_report = (
        adata.obs.groupby(leiden_col)[label_col]
        .first()
        .reset_index()
        .rename(columns={leiden_col: "cluster"})
    )

    output_excel_file = output_dir / "cluster_to_celltype_mapping.xlsx"
    df_report.to_excel(output_excel_file, index=False)

    sc.pl.umap(
        adata,
        color=label_col,
        save=f"_{leiden_col}_celltypeai_results.png",
    )

    multi_umap(adata, leiden_col, label_col, output_dir)

    # =====================
    # Save AnnData with scores
    # =====================
    output_h5ad = output_dir / "celltypeai_annotated.h5ad"
    adata.write(output_h5ad)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Launch CellTypeAI analysis and generate plots."
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

    parser.add_argument("--species", default="Human",
                        help="Species name passed to CellTypeAI (e.g. Human, Mouse)"
    )

    parser.add_argument("--tissue", required=True,
                        help="Tissue passed to CellTypeAI. Must be one of: Adrenal, "
                             "Blood, Brain, Eye, Heart, Immune system, Intestine, "
                             "Kidney, Liver, Lung, Muscle, Pancreas, Placenta, "
                             "Spleen, Stomach, Thymus, Skin."
    )

    parser.add_argument("--model", default="phi4:14b",
                        help="Ollama model tag already pulled locally "
                             "(`ollama pull <model>`). Must match a model served "
                             "by `ollama serve`."
    )

    parser.add_argument("--num_genes", type=int, default=200,
                        help="Number of top marker genes per cluster passed to CellTypeAI"
    )

    parser.add_argument("--n_iterations", type=int, default=3,
                        help="Number of ensemble prompt iterations (mode annotation is kept)"
    )

    parser.add_argument("--verbose", action="store_true",
                        help="Verbose CellTypeAI prompting (also dumps the engineered "
                             "prompts to the working directory)"
    )

    args = parser.parse_args()

    main(
        raw_data_file=args.input,
        verdict_file=args.verdict,
        output_dir=args.output,
        species=args.species,
        tissue=args.tissue,
        model=args.model,
        num_genes=args.num_genes,
        n_iterations=args.n_iterations,
        verbose=args.verbose,
    )
