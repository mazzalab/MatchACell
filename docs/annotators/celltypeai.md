# CellTypeAI

Annotates each Verdict cluster with [CellTypeAI](https://github.com/rhdaw/CellTypeAI),
which asks a local large language model, served by [Ollama](https://ollama.com),
to name the cell type from the cluster's marker genes. Nothing is sent to a
cloud service.

| | |
| --- | --- |
| Runs | when `celltypeai.tissue` is not empty (**on by default**) |
| Labels | per cluster |
| Rule / script | `workflow/rules/celltypeai.smk` → `workflow/scripts/celltypeai_tool.py` |
| Environment | `workflow/envs/celltypeai.yaml` (Python 3.12, `celltypeai` from PyPI, `ollama` from conda-forge) |
| Needs | nothing extra: the rule starts its own Ollama server and pulls the model if missing |

## Configuration

```yaml
matchacell_annotation:
  celltypeai:
    species: "Human"
    tissue: "Immune system"   # empty: skip CellTypeAI
    model: "phi3:mini"        # Ollama model tag; pulled into ollama_models if missing
    num_genes: 200            # marker genes per cluster in the prompt
    n_iterations: 1           # prompts per cluster; the most frequent answer is kept
    verbose: false            # true also writes the prompts to the working directory
    ollama_models: "/data/data_storage/research/Multiomics/ollama_models"  # empty: ~/.ollama/models
```

`tissue` must be one of: Adrenal, Blood, Brain, Eye, Heart, Immune system,
Intestine, Kidney, Liver, Lung, Muscle, Pancreas, Placenta, Spleen, Stomach,
Thymus, Skin.

**It's on by default.** The shipped `config.yaml` sets
`tissue: "Immune system"`, so `annotation` runs CellTypeAI. Set `tissue: ""`
to skip it.

## Ollama server

The rule's environment includes the `ollama` binary, and the rule manages the
server itself:

1. Starts `ollama serve` on a free port of the compute node, with
   `OLLAMA_MODELS` set to `ollama_models`. The server log goes to
   `annotation/CellTypeAI/ollama_serve.log`.
2. Waits up to 2 minutes for the server, and fails if it doesn't come up.
3. Runs `ollama pull <model>` if the model isn't in the store yet. This needs
   internet access from the node.
4. Points the script at the server through `OLLAMA_HOST`, then stops the
   server when the job ends.

CellTypeAI hardcodes `http://localhost:11434`. The script redirects those
requests to `OLLAMA_HOST`, so jobs on the same node don't collide.

If no cluster gets a label (for example because every request failed), the
script exits with an error instead of writing an all-NaN annotation.

Models trade accuracy for memory. `phi3:mini` needs about 2–3 GB of RAM;
`phi4:14b` about 9 GB; `qwen3:32b` or `qwen3:235b` are more accurate but much
heavier. Set `cluster_resources.celltypeai.mem_mb` to fit the model plus the
AnnData.

## How it works

1. Copies the Verdict clustering into a column named `leiden`, because
   CellTypeAI only looks for `leiden` or `louvain`. An existing `leiden` column
   is overwritten.
2. CellTypeAI ranks marker genes per cluster with a Wilcoxon test on the
   log-normalized `.X`. The script leaves `.raw` unset so that the ranking
   doesn't fall back to raw counts.
3. For each cluster, it prompts the model `n_iterations` times with the top
   `num_genes` markers, the species and the tissue, and keeps the most frequent
   answer.
4. The answer is stored on each cell as `cell_type_ai`.

## Outputs

In `annotation/CellTypeAI/`:

- `celltypeai_annotated.h5ad`: `cell_type_ai` per cell, plus the `leiden` copy.
- `cluster_to_celltype_mapping.xlsx`: `cell_type_ai` per cluster.
- UMAPs of `cell_type_ai` (overall and one panel per label).

## Notes

- Small models answer quickly but less reliably. Raising `n_iterations` smooths
  out their variability at the cost of runtime.
- The script's own defaults (`phi4:14b`, 3 iterations) only apply when it's run
  outside the workflow.
