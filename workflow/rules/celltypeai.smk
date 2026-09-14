# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: CellTypeAI
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config=_MC["celltypeai"]

rule celltypeai:
    """Run CellTypeAI annotation (local Ollama LLM ensemble)"""
    input:
        # Map each sample ID to its input .h5ad path from the config.
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"),
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CellTypeAI", "celltypeai_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","CellTypeAI"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        species=method_config["species"],
        tissue=method_config["tissue"],
        model=method_config["model"],
        num_genes=method_config["num_genes"],
        n_iterations=method_config["n_iterations"],
        # Precomputed here (plain python) because the shell: block below only
        # does simple {attribute} substitution, not inline python ternaries.
        verbose_flag=("--verbose" if method_config.get("verbose", False) else ""),
    threads: 8
    conda:
        "../envs/celltypeai.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""
        python workflow/scripts/celltypeai_tool.py \
            --input {input.h5ad} \
            --verdict {params.verdict_file} \
            --output {params.outdir} \
            --species {params.species:q} \
            --tissue {params.tissue:q} \
            --model {params.model:q} \
            --num_genes {params.num_genes} \
            --n_iterations {params.n_iterations} \
            {params.verbose_flag}
        """
