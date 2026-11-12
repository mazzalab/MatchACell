# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: CellWhisperer
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.
#
# CellWhisperer ships a pixi environment (PyTorch, Geneformer, its own forks),
# not a conda one, so this rule has no `conda:` directive: it runs the script
# with `pixi run` inside the CellWhisperer checkout given by `repo_dir`.

import os
import shlex

_CW =config["matchacell_annotation"].get("cellwhisperer", {})
_CW_REPO = _CW.get("repo_dir", "")
_CW_MODEL = _CW.get("model") or (
    os.path.join(_CW_REPO, "results", "models", "jointemb", "cellwhisperer_clip_v1.ckpt")
    if _CW_REPO else ""
)


rule cellwhisperer:
    """Annotate one clustered sample with CellWhisperer text-transcriptome scores."""
    input:
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"),
        verdict=os.path.join(outputDir, "results", "{sample}", "matchacell", "MatchA_Verdict.txt"),
        model=[_CW_MODEL] if _CW_REPO else [],
        manifest=[os.path.join(_CW_REPO, "pixi.toml")] if _CW_REPO else [],
        signatures=[] if _CW.get("labels") else [config["matchacell_annotation"]["annot_file"]],
    output:
        annotated_h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "annotation", "CellWhisperer", "cellwhisperer_annotated.h5ad"),
        mapping=os.path.join(outputDir, "results", "{sample}", "matchacell", "annotation", "CellWhisperer", "cluster_to_celltype_mapping.xlsx"),
    params:
        outdir=os.path.join(outputDir, "results", "{sample}", "matchacell", "annotation", "CellWhisperer"),
        pixi=_CW.get("pixi", "pixi"),
        repo_dir=_CW_REPO,
        # Functions, so Snakemake doesn't read "{label}" as a wildcard.
        labels=lambda wc: (
            "--labels " + " ".join(shlex.quote(str(x)) for x in _CW["labels"])
            if _CW.get("labels") else ""
        ),
        annot_file=config["matchacell_annotation"].get("annot_file", ""),
        template=lambda wc: _CW.get("template", "{label}"),
        threshold=_CW.get("threshold", 0.0),
        batch_size=_CW.get("batch_size", 32),
        device=_CW.get("device", "auto"),
        script=os.path.join(workflow.basedir, "workflow", "scripts", "cellwhisperer_tool.py"),
    threads: _CW.get("n_threads", 8)
    resources:
        mem_mb=cluster_resource("cellwhisperer", "mem_mb", 32000),
        runtime=cluster_resource("cellwhisperer", "runtime", 720),
    message:
        "MatchACell · CellWhisperer annotation · sample {wildcards.sample}"
    shell:
        r"""
        {params.pixi:q} run --manifest-path {input.manifest:q} python {params.script:q} \
            --input {input.h5ad:q} \
            --verdict {input.verdict:q} \
            --output {params.outdir:q} \
            --model {input.model:q} \
            --repo_dir {params.repo_dir:q} \
            --annot_file {params.annot_file:q} \
            --template {params.template:q} \
            --threshold {params.threshold} \
            --batch_size {params.batch_size} \
            --device {params.device:q} \
            --n_threads {threads} {params.labels}
        """
