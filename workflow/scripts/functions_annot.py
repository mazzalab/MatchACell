def rearrange_sign(raw_data_filtered, annot_file):

    import pandas as pd

    signature_dict = {
        sheet_name: [
            gene.strip()
            for gene in df.iloc[:, 0].dropna().astype(str)
            if gene.strip() in raw_data_filtered.var_names
        ]
        for sheet_name, df in pd.read_excel(
            annot_file, sheet_name=None, header=None
        ).items()
    }

    return signature_dict


def extract_best_res(verdict_file):
    
    import re

    with open(verdict_file, "r") as f:
        text = f.read()

    match = re.search(r"Recommended resolution\s*:\s*(\S+)", text)

    if match:
        return match.group(1)

    raise ValueError(
        f"'Recommended resolution' not found in {verdict_file}"
    )

def multi_umap(adata, col_base, col_name, output_umap):
    
    import matplotlib.pyplot as plt
    import scanpy as sc
    import numpy as np
    
    col_names = sorted(adata.obs[col_name].unique().tolist())
    n_types = len(col_names)

    # Total subplots: 1 (all samples) + 1 per sample
    ncols = 4
    nrows = int(np.ceil((n_types + 1) / ncols))
    _, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows), dpi=150)
    axes = axes.flatten()

    # --- First plot: All samples colored ---
    ax = axes[0]
    sc.pl.umap(adata, color=col_base, ax=ax, show=False, legend_loc=None, title="All Samples")
    
    # --- One plot per sample ---
    for i, sample in enumerate(col_names):
        ax = axes[i + 1]  # +1 because 0 is used
        mask = adata.obs[col_name] == sample

        # Background in gray
        ax.scatter(
            adata.obsm["X_umap"][~mask, 0],
            adata.obsm["X_umap"][~mask, 1],
            c="lightgrey",
            s=1,
            alpha=0.2,
            linewidths=0
        )

        # Foreground in color
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=f"C{i % 10}",
            s=1,
            alpha=0.8,
            linewidths=0
        )

        ax.set_title(sample, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_frame_on(False)

    # Hide unused axes
    for j in range(n_types + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.5, hspace=0.4)
    name_col=col_name.replace(" ","_")
    plt.savefig(output_umap / f"umap_per_{name_col}.png", dpi=150, bbox_inches="tight")