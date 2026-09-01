# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 0 — stage every sample as .h5ad
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.
#
# Every sample must land as results/{sample}/h5/{sample}.h5ad before Step 1
# (matchacell_cluster_stability) runs, regardless of what format it was
# supplied in:
#   - config["dtype"] == "rds"  -> workflow/scripts/rds2h5.R converts the
#     Seurat/SingleCellExperiment .rds into .h5ad (see the script's header
#     for why it writes h5ad directly via rhdf5 instead of SeuratDisk).
#   - anything else (native .h5ad) -> copied into place unchanged.
# The rest of the pipeline always reads from that h5/ location, so it never
# has to care which format the raw sample arrived in.

_DTYPE = config.get("dtype", "h5ad")


rule stage_h5:
    """Stage one sample as .h5ad under results/{sample}/h5/."""
    input:
        raw=lambda wc: config["samples"][wc.sample],
    output:
        h5ad=outputDir + "results/{sample}/h5/{sample}.h5ad",
    params:
        outdir=outputDir + "results/{sample}",
        dtype=_DTYPE,
    conda:
        # rds2h5.R needs Seurat/rhdf5/SingleCellExperiment/optparse -- the
        # same R env addmodulescore already uses, reused here rather than
        # duplicated into its own env file.
        "../envs/addmodulescore.yaml"
    message:
        "MatchACell · stage-h5 · sample {wildcards.sample}"
    shell:
        r"""
        if [ "{params.dtype}" = "rds" ]; then
            Rscript workflow/scripts/rds2h5.R \
                --input {input.raw} \
                --output {params.outdir}
            mv {params.outdir}/h5/*.h5ad {output.h5ad}
        else
            mkdir -p "$(dirname {output.h5ad})"
            cp {input.raw} {output.h5ad}
        fi
        """
