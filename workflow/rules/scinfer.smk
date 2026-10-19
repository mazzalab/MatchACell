# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: ScInfeR (reference mapping)
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir` and `config` are provided by the top-level Snakefile that
# includes this module.

import os

_MC = config["matchacell_annotation"]
method_config = _MC.get("scinfer", {})

rule scinfer:
    """Run ScInfeR reference-mapping annotation"""
    input:
        # Seurat object straight from Step 1's .rds companion (rule
        # matchacell_cluster_stability_rds), same as addmodulescore.
        rds=os.path.join(outputDir, "results", "{sample}", "matchacell", "clustered_multi_resolution.rds"),
        # ScInfeR.R needs a Seurat .rds; a .h5ad reference_file is converted
        # via h52rds.R first (below), same converter rule stage_h5 runs.
        reference=method_config.get("reference_file", ""),
    output:
        annotated_h5ad=os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "ScInfeR", "scinfer_annotated.h5ad")
    params:
        outdir=os.path.join(outputDir,"results","{sample}","matchacell","annotation","ScInfeR"),
        verdict_file=os.path.join(outputDir,"results","{sample}","matchacell","MatchA_Verdict.txt"),
        ref_label_col=method_config.get("ref_label_col", "cell_type"),
        tissue=method_config.get("tissue", ""),
        num_marker_per_ct=method_config.get("num_marker_per_ct", 10),
        own_weightage=method_config.get("own_weightage", 0.5),
        n_neighbor=method_config.get("n_neighbor", 10),
    threads: 4
    resources:
        mem_mb=cluster_resource("scinfer", "mem_mb", 32000),
        runtime=cluster_resource("scinfer", "runtime", 240),
    conda:
        "../envs/scinfer.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""
        REF_RDS={input.reference}
        if [[ "{input.reference}" == *.h5ad ]]; then
            Rscript workflow/scripts/h52rds.R \
                --input {input.reference} \
                --output {params.outdir}
            REF_BASENAME=$(basename {input.reference} .h5ad)
            REF_RDS={params.outdir}/rds/${{REF_BASENAME}}.rds
        fi

        Rscript workflow/scripts/ScInfeR.R \
            --input {input.rds} \
            --reference $REF_RDS \
            --ref_label_col {params.ref_label_col:q} \
            --tissue "{params.tissue}" \
            --num_marker_per_ct {params.num_marker_per_ct} \
            --own_weightage {params.own_weightage} \
            --n_neighbor {params.n_neighbor} \
            --verdict {params.verdict_file} \
            --output {params.outdir}

        Rscript workflow/scripts/rds2h5.R \
            --input {params.outdir}/rds/scinfer_annotated.rds \
            --output {params.outdir}
        mv {params.outdir}/h5/*.h5ad {output.annotated_h5ad}
        """
