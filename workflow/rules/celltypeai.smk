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
        # Shared Ollama model store (OLLAMA_MODELS); empty -> ~/.ollama/models.
        ollama_models=method_config.get("ollama_models", ""),
        # Thinking step of reasoning models (qwen3...): "false" skips it
        # (much faster), "true" keeps it, empty -> model default.
        think_flag=(
            "" if method_config.get("think") is None
            else "--think " + str(bool(method_config["think"])).lower()
        ),
    threads: 8
    resources:
        mem_mb=cluster_resource("celltypeai", "mem_mb", 16000),
        runtime=cluster_resource("celltypeai", "runtime", 480),
    conda:
        "../envs/celltypeai.yaml"
    message:
        "MatchACell · cluster-annotation · sample {wildcards.sample}"
    shell:
        r"""
        # Start a private Ollama server on a free port of this node, so the
        # rule doesn't depend on one already running (and two jobs on the
        # same node don't collide). celltypeai_tool.py reads OLLAMA_HOST.
        mkdir -p {params.outdir}
        if [ -n {params.ollama_models:q} ]; then
            export OLLAMA_MODELS={params.ollama_models:q}
        fi
        OLLAMA_PORT=$(python -c 'import socket; s = socket.socket(); s.bind(("127.0.0.1", 0)); print(s.getsockname()[1])')
        export OLLAMA_HOST=127.0.0.1:$OLLAMA_PORT
        OLLAMA_LOG={params.outdir}/ollama_serve.log
        ollama serve > "$OLLAMA_LOG" 2>&1 &
        OLLAMA_PID=$!
        trap 'kill $OLLAMA_PID 2>/dev/null || true' EXIT

        for i in $(seq 1 60); do
            if ollama list > /dev/null 2>&1; then break; fi
            sleep 2
        done
        if ! ollama list > /dev/null 2>&1; then
            echo "[CellTypeAI] ERROR: ollama serve did not start, see $OLLAMA_LOG" >&2
            exit 1
        fi
        echo "[CellTypeAI] ollama serve up on $OLLAMA_HOST (models: ${{OLLAMA_MODELS:-~/.ollama/models}})"

        # Serialize pulls into the shared store: two jobs pulling the same
        # model at once write the same -partial blob.
        mkdir -p "${{OLLAMA_MODELS:-$HOME/.ollama/models}}"
        (
            flock 9
            if ! ollama show {params.model:q} > /dev/null 2>&1; then
                echo "[CellTypeAI] Model {params.model} not found locally, pulling..."
                ollama pull {params.model:q}
            fi
        ) 9> "${{OLLAMA_MODELS:-$HOME/.ollama/models}}/.pull.lock"

        python workflow/scripts/celltypeai_tool.py \
            --input {input.h5ad} \
            --verdict {params.verdict_file} \
            --output {params.outdir} \
            --species {params.species:q} \
            --tissue {params.tissue:q} \
            --model {params.model:q} \
            --num_genes {params.num_genes} \
            --n_iterations {params.n_iterations} \
            --threads {threads} {params.think_flag} {params.verbose_flag}
        """
        # Optional flags share the last line (see cluster_stability.smk).
