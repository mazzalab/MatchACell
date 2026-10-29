#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

# --- Configurazione Backend Matplotlib (No Pop-up) ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import celltypeai as cta
import pandas as pd

import scanpy as sc

# Import dei moduli custom del workflow
from functions_annot import extract_best_res, multi_umap


def patch_ollama_requests(num_thread=None, think=None):
    """Adjust CellTypeAI's requests to Ollama.

    celltypeai hardcodes http://localhost:11434 and a fixed request body
    (prompt_submit.py), so patch requests.post to:
    - send them to the server in $OLLAMA_HOST, so the rule can run its own
      `ollama serve` on a free port;
    - set options.num_thread, otherwise llama-server uses every core of the
      node, not the threads PBS allotted;
    - set "think", so reasoning models (e.g. qwen3) can skip the thinking
      tokens that celltypeai discards anyway.
    """
    host = os.environ.get("OLLAMA_HOST", "").strip()
    if host:
        if not host.startswith(("http://", "https://")):
            host = "http://" + host
        host = host.rstrip("/")

    import celltypeai.prompt_submit as ps
    _post = ps.requests.post

    def _patched_post(url, *args, **kwargs):
        if host:
            url = url.replace("http://localhost:11434", host)
        payload = kwargs.get("json")
        if isinstance(payload, dict):
            if num_thread:
                payload.setdefault("options", {})["num_thread"] = num_thread
            if think is not None:
                payload["think"] = think
        return _post(url, *args, **kwargs)

    ps.requests.post = _patched_post
    print(f"[CellTypeAI] Ollama server: {host or 'http://localhost:11434'}, "
          f"num_thread: {num_thread or 'auto'}, think: {think}")


def main(raw_data_file, verdict_file, output_dir, species, tissue, model, num_genes, n_iterations, verbose,
         threads=None, think=None):

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
    #
    # celltypeai also names each rank_genes_groups group by its position
    # (str(0), str(1), ...) instead of its real label, so any label that is
    # not 0..n-1 (e.g. sub-clusters "15_0") gets the genes of another group
    # or ends up NaN. Feed it positional codes in natural order instead; the
    # per-cell cell_type_ai labels stay valid for the real leiden_col.
    labels = adata.obs[leiden_col].astype(str)
    ordered = sorted(
        labels.unique(),
        key=lambda x: [(0, int(p)) if p.isdigit() else (1, p) for p in x.split("_")]
    )
    codes = {label: str(i) for i, label in enumerate(ordered)}
    adata.obs["leiden"] = pd.Categorical(
        labels.map(codes), categories=[str(i) for i in range(len(ordered))]
    )

    patch_ollama_requests(num_thread=threads, think=think)

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

    # celltypeai swallows connection errors and leaves every label NaN;
    # fail here so Snakemake doesn't mark an empty annotation as done.
    if label_col not in adata.obs or adata.obs[label_col].isna().all():
        sys.exit("[CellTypeAI] ERROR: no cluster was annotated (all NaN). "
                 "Check that the Ollama server was reachable and the model "
                 "is available.")

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

    parser.add_argument("--threads", type=int, default=None,
                        help="CPU threads for the Ollama model (options.num_thread). "
                             "Default: Ollama decides (all cores of the node)."
    )

    parser.add_argument("--think", choices=["true", "false"], default=None,
                        help="Enable/disable the thinking step of reasoning models "
                             "(e.g. qwen3). Default: the model's own default."
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
        threads=args.threads,
        think=None if args.think is None else args.think == "true",
    )
