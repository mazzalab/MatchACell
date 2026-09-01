# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: AddModuleScore
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config=_MC["addmodulescore"]

rule addmodulescore:
    """Run AddModuleScore (Seurat) annotation"""
    input:
        # Seurat object straight from Step 1's .rds companion (rule
        # matchacell_cluster_stability_rds) -- no h5ad->Seurat conversion
        # needed here anymore, h52rds.R already did it once upstream.
        rds=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.rds"),
        signatures=_MC["annot_file"]
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "AddModuleScore", "addmodulescore_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","AddModuleScore"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        threshold=method_config["threshold"]
    threads: 4
    conda:
        "../envs/addmodulescore.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""
        Rscript workflow/scripts/addmodulescore_tool.R \
            --input {input.rds} \
            --annot {input.signatures} \
            --thr {params.threshold} \
            --verdict {params.verdict_file} \
            --output {params.outdir}

        Rscript workflow/scripts/rds2h5.R \
            --input {params.outdir}/rds/addmodulescore_annotated.rds \
            --output {params.outdir}
        mv {params.outdir}/h5/*.h5ad {output.annotated_h5ad}
        """
