# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 3 — cross-annotator consensus on the Cell Ontology
# ════════════════════════════════════════════════════════════════════════
#
# `outputDir`, `config` and the annotator switches (_CYTETYPE_TOKEN,
# _CELLTYPEAI_TISSUE, _SCANVI_REFERENCE, _SCPARADISE_MODEL) are provided by
# the top-level Snakefile that includes this module.

import shlex

_CONS = config.get("matchacell_consensus") or {}

# Pinned Cell Ontology release, downloaded once per output directory unless
# `ontology` points at a local copy (e.g. for compute nodes without internet).
CL_RELEASE = "2026-06-08"
_CL_OBO = _CONS.get("ontology") or outputDir + f"resources/cl-basic_{CL_RELEASE}.obo"

# Annotator name (as consensus_tool.py knows it) -> output file, for the
# annotators this config enables.
_ANNOTATOR_OUTPUTS = {
    "ScoreGenes": ("ScoreGenes", "score_genes", True),
    "AddModuleScore": ("AddModuleScore", "addmodulescore", True),
    "CIA": ("CIA", "cia", True),
    "CellTypist": ("CellTypist", "celltypist", True),
    "CyteType": ("CyteType", "cytetype", bool(_CYTETYPE_TOKEN)),
    "CellTypeAI": ("CellTypeAI", "celltypeai", bool(_CELLTYPEAI_TISSUE)),
    "scANVI": ("scANVI", "scanvi", bool(_SCANVI_REFERENCE)),
    "scParadise": ("scParadise", "scparadise", bool(_SCPARADISE_MODEL)),
}


def consensus_annotations(sample):
    exclude = set(_CONS.get("exclude") or [])
    return {
        name: os.path.join(outputDir, "results", sample, "matchacell", "annotation",
                           folder, f"{stem}_annotated.h5ad")
        for name, (folder, stem, enabled) in _ANNOTATOR_OUTPUTS.items()
        if enabled and name not in exclude
    }


rule fetch_cell_ontology:
    """Download the pinned Cell Ontology release (basic edition: is_a only)."""
    output:
        outputDir + "resources/cl-basic_{release}.obo",
    message:
        "MatchACell · consensus · Cell Ontology {wildcards.release}"
    shell:
        r"""
        curl -fsSL -o {output:q}.part \
            http://purl.obolibrary.org/obo/cl/releases/{wildcards.release}/cl-basic.obo
        mv {output:q}.part {output:q}
        """


rule matchacell_consensus:
    """Cross-annotator consensus for one sample at the Verdict resolution."""
    input:
        h5ad=outputDir + "results/{sample}/matchacell/clustered_multi_resolution.h5ad",
        verdict=outputDir + "results/{sample}/matchacell/MatchA_Verdict.txt",
        ontology=_CL_OBO,
        label_maps=_CONS.get("label_maps") or [],
        annotations=lambda wc: list(consensus_annotations(wc.sample).values()),
    output:
        table=outputDir + "results/{sample}/matchacell/consensus/consensus_by_cluster.csv",
        h5ad=outputDir + "results/{sample}/matchacell/consensus/consensus_annotated.h5ad",
    params:
        outdir=outputDir + "results/{sample}/matchacell/consensus",
        annotations=lambda wc: " ".join(
            "-a " + shlex.quote(f"{name}={path}")
            for name, path in consensus_annotations(wc.sample).items()
        ),
        label_maps=" ".join("--label-map " + shlex.quote(p) for p in _CONS.get("label_maps") or []),
        min_agreement=_CONS.get("min_agreement", 0.5),
    resources:
        mem_mb=cluster_resource("matchacell_consensus", "mem_mb", 16000),
        runtime=cluster_resource("matchacell_consensus", "runtime", 60),
    conda:
        "../envs/matchacell.yaml"
    message:
        "MatchACell · consensus · sample {wildcards.sample}"
    shell:
        r"""
        python workflow/scripts/consensus_tool.py \
            --input {input.h5ad:q} \
            --verdict {input.verdict:q} \
            --ontology {input.ontology:q} \
            --min-agreement {params.min_agreement} \
            --output {params.outdir:q} {params.annotations} {params.label_maps}
        """
