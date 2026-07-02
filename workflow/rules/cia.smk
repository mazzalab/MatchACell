# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: CIA
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config=_MC["cia"]

rule cia:
    """Run CIA annotation"""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"), #lambda wc: config["samples"][wc.sample],
        signatures=_MC["annot_file"]
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CIA", "cia_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","CIA"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        threshold=method_config["threshold"],
        ncpus=method_config["ncpus"]
    threads: 8
    conda:
        "../envs/cia.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""        
        python workflow/scripts/cia_tool.py \
            --input {input.h5ad} \
            --annot {input.signatures} \
            --verdict {params.verdict_file} \
            --output {params.outdir} \
            --thr {params.threshold} \
            --ncpus {params.ncpus}
        """
