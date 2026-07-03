# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: CellTypist
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config=_MC["celltypist"]

rule celltypist:
    """Run CIA annotation"""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"), #lambda wc: config["samples"][wc.sample],
        signatures=_MC["annot_file"]
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CellTypist", "celltypist_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","CellTypist"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        models=method_config["models"],
    threads: 8
    conda:
        "../envs/celltypist.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""        
        python workflow/scripts/celltypist_tool.py \
            --input {input.h5ad} \
            --verdict {params.verdict_file} \
            --models_list {params.models} \
            --output {params.outdir}
        """
