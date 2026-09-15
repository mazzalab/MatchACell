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
| Needs | a running Ollama server with the configured model pulled |

## Configuration

```yaml
matchacell_annotation:
  celltypeai:
    species: "Human"
    tissue: "Immune system"   # empty: skip CellTypeAI
    model: "phi3:mini"        # an Ollama model that is already pulled
    num_genes: 200            # marker genes per cluster in the prompt
    n_iterations: 1           # prompts per cluster; the most frequent answer is kept
    verbose: false            # true also writes the prompts to the working directory
```

`tissue` must be one of: Adrenal, Blood, Brain, Eye, Heart, Immune system,
Intestine, Kidney, Liver, Lung, Muscle, Pancreas, Placenta, Spleen, Stomach,
Thymus, Skin.

**It's on by default.** The shipped `config.yaml` sets
`tissue: "Immune system"`, so `annotation` runs CellTypeAI and fails if no
Ollama server is available. Set `tissue: ""` to skip it.

## Setting up Ollama

The rule's environment includes the `ollama` binary. Start the server and pull
the model where the rule will run:

```bash
ollama serve &          # keep it running
ollama pull phi3:mini
```

Models trade accuracy for memory. `phi3:mini` needs about 2–3 GB of RAM;
`phi4:14b`, `qwen3:32b` or `qwen3:235b` are more accurate but much heavier. On a
compute cluster, the server must be reachable from the node that runs the rule.

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
