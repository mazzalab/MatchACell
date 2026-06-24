import os

# ════════════════════════════════════════════════════════════════════════
#  MatchACell — Snakemake entrypoint
#  Brewed with MatchACell · matcha-grade single-cell consensus
# ════════════════════════════════════════════════════════════════════════

# Normalize the output directory so concatenation always has a separator.
outputDir = config["output_dir"].rstrip("/") + "/"

# ── Included modules ──────────────────────────────────────────────────────
# Step 1. Future annotator steps will be added as additional includes here.
include: "workflow/rules/cluster_stability.smk"


SAMPLES = (
    list(config["samples"].keys())
    if isinstance(config["samples"], dict)
    else list(config["samples"])
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
