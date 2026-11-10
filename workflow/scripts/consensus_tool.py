#!/usr/bin/env python3
"""
MatchACell - Step 3: cross-annotator consensus on the Cell Ontology.

For the clusters of the MatchA Verdict resolution, every annotator votes with
its per-cell labels: in each cluster it gives each label the fraction of the
cluster's cells carrying it (Unknown/Unassigned cells abstain, ties such as
"CD8 T, CD4 T" split the cell between the labels). Labels are mapped to Cell
Ontology (CL) terms, and a vote for a term also supports every is_a ancestor of
it: "CD4-positive, alpha-beta T cell" and "CD8-positive, alpha-beta T cell"
disagree, but both support "T cell".

The consensus of a cluster is the most specific term supported by more than
`--min-agreement` of the annotators that voted on it:
  agreed   the consensus term is a label some annotator gave;
  partial  annotators only agree on a coarser common ancestor;
  none     no term reaches the bar (e.g. unmapped free-text labels).
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import anndata as ad
import pandas as pd

from cell_ontology import CellOntology, normalize
from functions_annot import extract_best_res

# Per-cell label column(s) of each annotator. CellTypist writes one column per
# model; its models share the annotator's single vote.
LABEL_COLUMNS = {
    "ScoreGenes": ["cell_type_pred"],
    "AddModuleScore": ["cell_type_pred"],
    "CIA": ["CIA prediction default"],
    "CellTypist": "celltypist_label_majority_voting_",
    "scANVI": ["cell_type_pred"],
    "scParadise": ["cell_type_pred"],
    "CyteType": ["cytetype_annotation_{leiden}"],
    "CellTypeAI": ["cell_type_ai"],
}
# CL ids the annotator gives itself, used before any label mapping.
ID_COLUMNS = {"CyteType": "cytetype_cellOntologyTermID_{leiden}"}
# Annotators that join tied labels with ", " (marker-score ties).
TIED = {"ScoreGenes", "AddModuleScore"}
ABSTAIN = {"", "unknown", "unassigned", "nan", "none", "na", "n/a", "unclassified"}
_CL_ID = re.compile(r"CL:\d{7}")


# =========================================================================== #
# Labels -> Cell Ontology terms
# =========================================================================== #
class LabelMapper:
    """label -> (term, how). Label maps are TSVs with `label` and `cl_id`
    columns; a later map overrides an earlier one. Unmapped labels become a
    free-text term "free:<label>" with no ancestors."""

    def __init__(self, cl: CellOntology, label_maps=()):
        self.cl, self.table = cl, {}
        for path in label_maps:
            df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
            unknown = sorted(set(df["cl_id"]) - set(cl.names))
            if unknown:
                raise ValueError(f"{path}: CL ids not in the ontology: {unknown}")
            self.table.update({normalize(lab): cid for lab, cid in zip(df["label"], df["cl_id"])})

    def __call__(self, label, cl_id=None):
        if cl_id and _CL_ID.fullmatch(str(cl_id).strip()) and str(cl_id).strip() in self.cl:
            return str(cl_id).strip(), "annotator"
        key = normalize(label)
        if key in self.table:
            return self.table[key], "label map"
        if _CL_ID.fullmatch(str(label).strip()) and str(label).strip() in self.cl:
            return str(label).strip(), "CL id"
        term = self.cl.lookup(label)
        if term:
            return term, "CL name"
        return f"free:{key}", "unmapped"

    def name(self, term):
        return self.cl.names.get(term, term.removeprefix("free:"))

    def ancestors(self, term):
        return self.cl.ancestors(term) if term in self.cl else frozenset([term])

    def specificity(self, term):
        return self.cl.specificity(term) if term in self.cl else 0


# =========================================================================== #
# Votes
# =========================================================================== #
def label_columns(method, obs, leiden_col):
    spec = LABEL_COLUMNS[method]
    if isinstance(spec, str):
        return [c for c in obs.columns if c.startswith(spec)]
    return [c.format(leiden=leiden_col) for c in spec if c.format(leiden=leiden_col) in obs]


def cluster_shares(labels: pd.Series, clusters: pd.Series, tied=False):
    """Long table (cluster, label, share): the fraction of each cluster's cells
    carrying each label. Abstaining cells count in the cluster size only."""
    df = pd.DataFrame({"cluster": clusters.astype(str).values,
                       "label": labels.astype(str).str.strip().values})
    size = df.groupby("cluster").size()
    if tied:
        df["label"] = df["label"].str.split(r"\s*,\s*")
        df["w"] = 1.0 / df["label"].str.len()
        df = df.explode("label")
    else:
        df["w"] = 1.0
    df = df[~df["label"].str.lower().isin(ABSTAIN)]
    out = df.groupby(["cluster", "label"], as_index=False)["w"].sum()
    out["share"] = out["w"] / out["cluster"].map(size)
    return out.drop(columns="w")


def collect_votes(annotations: dict, clusters: pd.Series, leiden_col, mapper: LabelMapper):
    """Long vote table: annotator, voter, cluster, label, cl_id, mapped_by,
    share, weight (share x the voter's part of the annotator's vote)."""
    frames = []
    for method, path in annotations.items():
        obs = ad.read_h5ad(path, backed="r").obs
        common = clusters.index.intersection(obs.index)
        if len(common) < 0.9 * len(clusters):
            raise ValueError(f"{method}: only {len(common)}/{len(clusters)} cells match "
                             f"the clustered object ({path})")
        cols = label_columns(method, obs, leiden_col)
        if not cols:
            print(f"[consensus] WARNING: {method} has no label column; skipped.")
            continue
        id_col = ID_COLUMNS.get(method, "").format(leiden=leiden_col)
        ids = (obs.loc[common, id_col].astype(str) if id_col in obs else None)
        for col in cols:
            shares = cluster_shares(obs.loc[common, col], clusters.loc[common],
                                    tied=method in TIED)
            if ids is not None:  # annotator-given CL id of each label
                first_id = pd.DataFrame({"label": obs.loc[common, col].astype(str).str.strip(),
                                         "id": ids}).groupby("label")["id"].first()
                shares["given_id"] = shares["label"].map(first_id)
            else:
                shares["given_id"] = None
            mapped = [mapper(lab, gid) for lab, gid in zip(shares["label"], shares["given_id"])]
            shares["cl_id"] = [m[0] for m in mapped]
            shares["mapped_by"] = [m[1] for m in mapped]
            shares["annotator"] = method
            shares["voter"] = (method if len(cols) == 1
                               else f"{method}:{col.removeprefix(LABEL_COLUMNS[method]).removesuffix('.pkl')}")
            shares["weight"] = shares["share"] / len(cols)
            frames.append(shares.drop(columns="given_id"))
    if not frames:
        raise ValueError("No annotator produced labels.")
    votes = pd.concat(frames, ignore_index=True)
    votes["cl_name"] = votes["cl_id"].map(mapper.name)
    return votes[["annotator", "voter", "cluster", "label", "cl_id", "cl_name",
                  "mapped_by", "share", "weight"]]


def consensus_for_cluster(votes: pd.DataFrame, mapper: LabelMapper, min_agreement=0.5):
    """Consensus of one cluster's votes (rows of `collect_votes`)."""
    n_annot = votes["annotator"].nunique()
    if n_annot == 0:
        return {"status": "none", "n_annotators": 0}
    leaf = votes.groupby("cl_id")["weight"].sum()
    support = defaultdict(float)
    for term, w in leaf.items():
        for a in mapper.ancestors(term):
            support[a] += w
    bar = min_agreement * n_annot
    passing = [t for t, s in support.items() if s > bar + 1e-12]
    top = leaf.idxmax()
    if passing:
        best = max(passing, key=lambda t: (mapper.specificity(t), support[t], t))
        status = "agreed" if best in leaf.index else "partial"
    else:
        best, status = top, "none"
    coarser = sum(w for t, w in leaf.items()
                  if t != best and t in mapper.ancestors(best))
    # Show the consensus in the annotators' own words when one of them used it.
    named = votes[votes["cl_id"] == best].groupby("label")["weight"].sum()
    return {
        "consensus_label": named.idxmax() if len(named) else mapper.name(best),
        "consensus_cl_id": best if best in mapper.cl else "",
        "consensus_cl_name": mapper.name(best),
        "agreement": support[best] / n_annot,
        # votes that don't contradict the consensus (it or a coarser term)
        "compatibility": (support[best] + coarser) / n_annot,
        "status": status,
        "top_label": mapper.name(top),
        "top_label_agreement": leaf[top] / n_annot,
        "n_annotators": n_annot,
    }


def consensus_table(votes, mapper, clusters: pd.Series, min_agreement=0.5, stability=None):
    rows = []
    order = sorted(clusters.astype(str).unique(),
                   key=lambda x: [(0, int(p)) if p.isdigit() else (1, p) for p in x.split("_")])
    size = clusters.astype(str).value_counts()
    for cl in order:
        v = votes[votes["cluster"] == cl]
        row = {"cluster": cl, "n_cells": int(size[cl])}
        if stability is not None:
            row["stability"] = stability.get(cl)
        row.update(consensus_for_cluster(v, mapper, min_agreement))
        # each voter's leading label, e.g. "CD4 T (0.93)"
        for voter, g in v.groupby("voter"):
            r = g.loc[g["share"].idxmax()]
            row[f"call:{voter}"] = f"{r['label']} ({r['share']:.2f})"
        rows.append(row)
    return pd.DataFrame(rows)


def read_stability(outdir: Path, leiden_col):
    res = leiden_col.replace("leiden_", "")
    csv = outdir / "stability" / "csv" / f"stability_{res}.csv"
    if csv.exists():
        df = pd.read_csv(csv, dtype={"Cluster": str})
    else:
        xlsx = outdir / "stability" / "Cluster_Stability_Summary.xlsx"
        if not xlsx.exists():
            return None
        df = pd.read_excel(xlsx, sheet_name=leiden_col[:31], dtype={"Cluster": str})
    return dict(zip(df["Cluster"], df["Jaccard_Stability"]))


# =========================================================================== #
# Main
# =========================================================================== #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("-i", "--input", required=True, help="Step 1 clustered .h5ad")
    p.add_argument("-v", "--verdict", required=True, help="MatchA_Verdict.txt")
    p.add_argument("-a", "--annotation", action="append", default=[], metavar="NAME=H5AD",
                   help=f"Annotator output; NAME in {', '.join(LABEL_COLUMNS)}.")
    p.add_argument("--ontology", required=True, help="Cell Ontology OBO (cl-basic.obo)")
    p.add_argument("--label-map", action="append", default=[],
                   help="TSV with label and cl_id columns; later maps take precedence.")
    p.add_argument("--min-agreement", type=float, default=0.5,
                   help="Share of the voting annotators a term needs more than.")
    p.add_argument("-o", "--output", required=True, help="Output directory")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    annotations = {}
    for item in args.annotation:
        name, _, path = item.partition("=")
        if name not in LABEL_COLUMNS:
            raise SystemExit(f"Unknown annotator {name!r}; expected one of {list(LABEL_COLUMNS)}")
        annotations[name] = path

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    leiden_col = extract_best_res(args.verdict)
    adata = ad.read_h5ad(args.input)
    clusters = adata.obs[leiden_col].astype(str)

    mapper = LabelMapper(CellOntology(args.ontology), args.label_map)
    votes = collect_votes(annotations, clusters, leiden_col, mapper)
    stability = read_stability(Path(args.input).parent, leiden_col)
    table = consensus_table(votes, mapper, clusters, args.min_agreement, stability)

    votes.to_csv(out / "annotation_long.csv", index=False)
    table.to_csv(out / "consensus_by_cluster.csv", index=False)
    cells = votes.assign(n_cells=votes["share"] * votes["cluster"].map(clusters.value_counts()))
    (cells.groupby(["annotator", "label", "cl_id", "cl_name", "mapped_by"], as_index=False)
          ["n_cells"].sum().round({"n_cells": 1})
          .sort_values(["mapped_by", "annotator", "n_cells"], ascending=[False, True, False])
          .to_csv(out / "label_vocabulary.csv", index=False))

    per_cluster = table.set_index("cluster")
    for col, key in [("consensus_label", "consensus_label"),
                     ("consensus_cl_id", "consensus_cl_id"),
                     ("consensus_cl_name", "consensus_cl_name"),
                     ("agreement", "consensus_agreement"),
                     ("status", "consensus_status")]:
        adata.obs[key] = clusters.map(per_cluster[col]).values
    adata.write_h5ad(out / "consensus_annotated.h5ad")

    n_unmapped = votes.loc[votes["mapped_by"] == "unmapped", "label"].nunique()
    print(f"[consensus] {leiden_col}: {len(table)} clusters, {len(annotations)} annotators; "
          f"status {table['status'].value_counts().to_dict()}; "
          f"{n_unmapped} unmapped label(s), see label_vocabulary.csv")


if __name__ == "__main__":
    main()
