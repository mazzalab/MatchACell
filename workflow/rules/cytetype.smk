# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: CyteType
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config=_MC["cytetype"]

rule cytetype:
    """Run CyteType annotation"""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"), #lambda wc: config["samples"][wc.sample],
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CyteType", "cytetype_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","CyteType"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        study_context=method_config["study_context"],
        title=method_config["title"],
        run_label=method_config["run_label"],
        n_top_genes=method_config["n_top_genes"],
        auth_token=method_config.get("api_token", ""),
    threads: 8
    conda:
        "../envs/cytetype.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""
        python workflow/scripts/cytetype_tool.py \
            --input {input.h5ad} \
            --verdict {params.verdict_file} \
            --output {params.outdir} \
            --study_context {params.study_context:q} \
            --title {params.title:q} \
            --run_label {params.run_label:q} \
            --n_top_genes {params.n_top_genes} \
            --auth_token {params.auth_token:q}
        """
