# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: scANVI (reference mapping)
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config = _MC.get("scanvi", {})

rule scanvi:
    """Run scANVI reference-mapping annotation"""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"), #lambda wc: config["samples"][wc.sample],
        reference=method_config.get("reference_file", ""),
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "scANVI", "scanvi_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","scANVI"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        label_col=method_config.get("label_col", "free_annotation"),
        unlabeled_category=method_config.get("unlabeled_category", "Unknown"),
        threshold=method_config.get("threshold", 0.7),
        max_epochs_ref=method_config.get("max_epochs_ref", 20),
        max_epochs_scanvi=method_config.get("max_epochs_scanvi", 20),
        max_epochs_query=method_config.get("max_epochs_query", 10),
        batch_size=method_config.get("batch_size", 2048),
        n_threads=method_config.get("n_threads", 10),
    threads: 10
    conda:
        "../envs/scanvi.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""
        python workflow/scripts/scanvi_tool.py \
            --input {input.h5ad} \
            --verdict {params.verdict_file} \
            --reference {input.reference} \
            --output {params.outdir} \
            --label_col {params.label_col:q} \
            --unlabeled_category {params.unlabeled_category:q} \
            --thr {params.threshold} \
            --max_epochs_ref {params.max_epochs_ref} \
            --max_epochs_scanvi {params.max_epochs_scanvi} \
            --max_epochs_query {params.max_epochs_query} \
            --batch_size {params.batch_size} \
            --n_threads {params.n_threads}
        """
