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
        # Always read the staged .h5ad (rule stage_h5 in rds2h5.smk), not
        # config["samples"] directly -- that path may point at a .rds.
        h5ad=outputDir + "results/{sample}/h5/{sample}.h5ad",
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


rule matchacell_cluster_stability_rds:
    """Companion .rds (Seurat object) of the clustered .h5ad, via h52rds.R --
    for downstream R-side use outside this pipeline."""
    input:
        h5ad=outputDir + "results/{sample}/matchacell/clustered_multi_resolution.h5ad",
    output:
        rds=outputDir + "results/{sample}/matchacell/clustered_multi_resolution.rds",
    params:
        outdir=outputDir + "results/{sample}/matchacell",
    conda:
        # h52rds.R needs zellkonverter (not just Seurat/rhdf5) -- same R env
        # addmodulescore already uses, reused here rather than duplicated.
        "../envs/addmodulescore.yaml"
    message:
        "MatchACell · cluster-stability (rds) · sample {wildcards.sample}"
    shell:
        r"""
        Rscript workflow/scripts/h52rds.R \
            --input {input.h5ad} \
            --output {params.outdir}
        mv {params.outdir}/rds/*.rds {output.rds}
        """
