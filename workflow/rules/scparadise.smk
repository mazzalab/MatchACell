# MatchACell · Step 2 — Cell Type Annotation: scParadise (scAdam)
import os

_SCP = config["matchacell_annotation"].get("scparadise", {})


rule scparadise:
    """Annotate one clustered sample with a local scAdam v2 model."""
    input:
        h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.h5ad"),
        verdict=os.path.join(outputDir, "results", "{sample}", "matchacell", "MatchA_Verdict.txt"),
        model=(
            [os.path.join(_SCP["model_dir"], "model_v2.pth")]
            if _SCP.get("model_dir") else []
        ),
        detector=(
            [os.path.join(_SCP["model_dir"], "unknown_detector.json")]
            if _SCP.get("model_dir") and _SCP.get("detect_unknown", False) else []
        ),
    output:
        annotated_h5ad=os.path.join(outputDir, "results", "{sample}", "matchacell", "annotation", "scParadise", "scparadise_annotated.h5ad"),
        mapping=os.path.join(outputDir, "results", "{sample}", "matchacell", "annotation", "scParadise", "cluster_to_celltype_mapping.xlsx"),
    params:
        outdir=os.path.join(outputDir, "results", "{sample}", "matchacell", "annotation", "scParadise"),
        model_dir=_SCP.get("model_dir", ""),
        layer=_SCP.get("layer", "lognorm"),
        label_level=_SCP.get("label_level", 0),
        threshold=_SCP.get("threshold", 0.0),
        unknown_flag="--detect_unknown" if _SCP.get("detect_unknown", False) else "",
        unknown_method=_SCP.get("unknown_method", "voting"),
        unknown_threshold=_SCP.get("unknown_threshold", 2.0),
        batch_size=_SCP.get("batch_size", 256),
        device=_SCP.get("device", "auto"),
        seed=_SCP.get("seed", 0),
    threads: _SCP.get("n_threads", 8)
    conda:
        "../envs/scparadise.yaml"
    message:
        "MatchACell · scParadise annotation · sample {wildcards.sample}"
    shell:
        r"""
        python workflow/scripts/scparadise_tool.py \
            --input {input.h5ad:q} \
            --verdict {input.verdict:q} \
            --model_dir {params.model_dir:q} \
            --output {params.outdir:q} \
            --layer={params.layer:q} \
            --label_level {params.label_level} \
            --thr {params.threshold} \
            {params.unknown_flag} \
            --unknown_method {params.unknown_method:q} \
            --unknown_threshold {params.unknown_threshold} \
            --batch_size {params.batch_size} \
            --device {params.device:q} \
            --n_threads {threads} \
            --seed {params.seed}
        """
