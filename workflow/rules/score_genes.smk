# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: ScoreGenes
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config=_MC["score_genes"]

rule score_genes:
    """Run ScoreGenes annotation"""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"), #lambda wc: config["samples"][wc.sample],
        signatures=_MC["annot_file"]
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "ScoreGenes", "score_genes_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","ScoreGenes"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        threshold=method_config["threshold"]
    threads: 8
    conda:
        "../envs/matchacell.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""        
        python workflow/scripts/score_genes.py \
            --input {input.h5ad} \
            --annot {input.signatures} \
            --thr {params.threshold} \
            --verdict {params.verdict_file} \
            --output {params.outdir}
        """
