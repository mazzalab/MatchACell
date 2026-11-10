#!/usr/bin/env python3
"""Build workflow/resources/label_maps/celltypist_immune.tsv.

Low-hierarchy labels of the Immune_All models take the Cell Ontology ID that
the CellTypist encyclopedia curates for them. High-hierarchy labels reuse the
ID of the low label with the same name, or else get the most specific CL term
that is an is_a ancestor of all their low labels.

    python tools/build_celltypist_label_map.py cl-basic.obo \
        workflow/resources/label_maps/celltypist_immune.tsv
"""
import argparse
import sys
from functools import reduce
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from cell_ontology import CellOntology  # noqa: E402

TABLE = ("https://raw.githubusercontent.com/Teichlab/celltypist_wiki/main/atlases/"
         "Pan_Immune_CellTypist/v2/tables/Basic_celltype_information.xlsx")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("obo")
    p.add_argument("output")
    p.add_argument("--table", default=TABLE, help="CellTypist encyclopedia table (URL or path).")
    args = p.parse_args()

    cl = CellOntology(args.obo)
    df = pd.read_excel(args.table)
    df["Cell Ontology ID"] = df["Cell Ontology ID"].astype(str).str.strip()
    missing = sorted(set(df["Cell Ontology ID"]) - set(cl.names))
    if missing:
        sys.exit(f"IDs not in this CL release (obsolete or newer): {missing}")

    rows = [(low, cid, "CellTypist encyclopedia")
            for low, cid in zip(df["Low-hierarchy cell types"], df["Cell Ontology ID"])]
    low_ids = dict(zip(df["Low-hierarchy cell types"], df["Cell Ontology ID"]))
    for high, g in df.groupby("High-hierarchy cell types"):
        if high in low_ids:
            continue
        common = reduce(frozenset.intersection, (cl.ancestors(i) for i in g["Cell Ontology ID"]))
        best = max(common, key=lambda t: (cl.specificity(t), t))
        rows.append((high, best, "common ancestor of its low-hierarchy labels"))

    out = pd.DataFrame(rows, columns=["label", "cl_id", "source"])
    out.insert(2, "cl_name", out["cl_id"].map(cl.names))
    out.to_csv(args.output, sep="\t", index=False)
    print(f"{len(out)} labels -> {args.output}")


if __name__ == "__main__":
    main()
