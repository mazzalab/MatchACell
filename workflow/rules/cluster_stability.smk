# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 1 — data-driven QC + Leiden cluster-stability optimizer
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

_MC = config["matchacell_cluster_stability"]


rule matchacell_cluster_stability:
    """Run QC + multi-resolution Leiden + bootstrap-Jaccard stability for one
    sample and emit the clustered .h5ad plus all diagnostic outputs."""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=lambda wc: config["samples"][wc.sample],
    output:
        clustered_h5ad=outputDir + "results/{sample}/matchacell/clustered_multi_resolution.h5ad",
    params:
        outdir=outputDir + "results/{sample}/matchacell",
        backend=_MC["backend"],
        n_iter=_MC["n_iter"],
        # Free-form extra flags forwarded verbatim to the script
        # (e.g. "--resolutions 0.2 0.5 1.0 --skip-tsne").
        extra=_MC.get("extra", ""),
    threads: 8
    conda:
        "../envs/matchacell.yaml"
    message:
        "MatchACell · cluster-stability · sample {wildcards.sample}"
    shell:
        r"""
        python workflow/scripts/matchacell_cluster_stability.py \
            --input   {input.h5ad} \
            --outdir  {params.outdir} \
            --backend {params.backend} \
            --n-iter  {params.n_iter} \
            {params.extra}
        """
