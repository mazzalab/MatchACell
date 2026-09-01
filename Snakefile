import os

# ════════════════════════════════════════════════════════════════════════
#  MatchACell — Snakemake entrypoint
#  Brewed with MatchACell · matcha-grade single-cell consensus
# ════════════════════════════════════════════════════════════════════════

# Normalize the output directory so concatenation always has a separator.
outputDir = config["output_dir"].rstrip("/") + "/"

# ── Included modules ──────────────────────────────────────────────────────
# Step 0. Stage every sample as .h5ad (converting from .rds when needed).
include: "workflow/rules/rds2h5.smk"

# Step 1. Future annotator steps will be added as additional includes here.
include: "workflow/rules/cluster_stability.smk"

# Step 2. Annotators.
include: "workflow/rules/score_genes.smk"
include: "workflow/rules/cia.smk"
include: "workflow/rules/celltypist.smk"
include: "workflow/rules/cytetype.smk"
include: "workflow/rules/addmodulescore.smk"


SAMPLES = (
    list(config["samples"].keys())
    if isinstance(config["samples"], dict)
    else list(config["samples"])
)

# CyteType needs an account/API token (see `cytetype setup`); skip it out of
# rule annotation's required inputs when none is configured, instead of
# failing the whole Step 2 run. Checked in config.yaml first (kept empty
# there, commit-safe) and then the CYTETYPE_API_TOKEN env var, so the token
# never has to be written into a tracked file.
_CYTETYPE_TOKEN = (
    config["matchacell_annotation"].get("cytetype", {}).get("api_token", "")
    or os.environ.get("CYTETYPE_API_TOKEN", "")
)


onstart:
    print("\033[1m\033[38;2;122;182;97m  Brewing MatchACell …\033[0m")


onsuccess:
    print("\033[1m\033[38;2;74;124;60m  MatchACell steeped successfully. 🍵\033[0m".replace(" 🍵", ""))


# ── Aggregating targets ───────────────────────────────────────────────────
rule cluster_stability:
    """Step 1 — QC + Leiden cluster stability for every sample."""
    input:
        expand(
            outputDir + "results/{sample}/matchacell/clustered_multi_resolution.h5ad",
            sample=SAMPLES,
        ),
        expand(
            outputDir + "results/{sample}/matchacell/clustered_multi_resolution.rds",
            sample=SAMPLES,
        ),

rule annotation:
    """Step 2 — Annotation."""
    input:
        score_genes=expand(
            os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "ScoreGenes", "score_genes_annotated.h5ad"),
            sample=SAMPLES,
        ),
        cia=expand(
            os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CIA", "cia_annotated.h5ad"),
            sample=SAMPLES,
        ),
        celltypist=expand(
            os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CellTypist", "celltypist_annotated.h5ad"),
            sample=SAMPLES,
        ),
        cytetype=(
            expand(
                os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "CyteType", "cytetype_annotated.h5ad"),
                sample=SAMPLES,
            )
            if _CYTETYPE_TOKEN else []
        ),
        addmodulescore=expand(
            os.path.join(outputDir, "results","{sample}", "matchacell", "annotation", "AddModuleScore", "addmodulescore_annotated.h5ad"),
            sample=SAMPLES,
        ),

    # input:
    #     lambda wildcards: [
    #         # Verifica in modo robusto che il file sia configurato e non sia una stringa vuota
    #         os.path.join(outputDir, "results", smp, "matchacell", "annotation", "ScoreGenes", "score_genes_annotated.h5ad")
    #         if (
    #             config.get("matchacell_annotation", {}).get("score_genes", {}).get("annot_file") 
    #             and str(config["matchacell_annotation"]["score_genes"]["annot_file"]).strip() != ""
    #         )
            
    #         # Altrimenti salta score_genes e usa lo step 1
    #         else os.path.join(outputDir, "results", smp, "matchacell", "clustered_multi_resolution.h5ad")
            
    #         for smp in SAMPLES
    #     ]
    
# rule all:
#     input:
#         rules.annotation.input