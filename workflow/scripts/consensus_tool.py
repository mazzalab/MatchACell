#!/usr/bin/env python3
"""
MatchACell - Step 3: cross-annotator consensus on the Cell Ontology.

For the clusters of the MatchA Verdict resolution, every annotator votes with
its per-cell labels, mapped to Cell Ontology (CL) terms: in each cluster a term
gets the fraction of the cluster's cells that carry it (Unknown/Unassigned cells
abstain). A vote for a term also supports every is_a ancestor of it:
"CD4-positive, alpha-beta T cell" and "CD8-positive, alpha-beta T cell"
disagree, but both support "T cell". This is the ontology vote of popV (Ergen
et al., Nat Genet 2024), applied to clusters.

The consensus of a cluster is the most specific term supported by more than
`--min-agreement` of the (weighted) annotators that voted on it:
  agreed   the consensus term is a label some annotator gave;
  partial  annotators only agree on a coarser common ancestor;
  none     no term reaches the bar (e.g. unmapped free-text labels).

With `--truth-col`, the annotators and the consensus are also scored per cell
against known labels in the Step 1 object.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from cell_ontology import CellOntology, normalize
from functions_annot import extract_best_res

# Per-cell label column(s) of each annotator. CellTypist writes one column per
# model; its models are combined into the annotator's single vote.
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
# Annotators whose label lists every signature above their threshold, best
# score first ("CD14 Mono, FCGR3A Mono"): their call is the first one.
RANKED = {"ScoreGenes", "AddModuleScore"}
# Sources of evidence. Annotators of one family share their inputs (the same
# marker workbook, reference atlases, or an LLM's knowledge of markers) and
# agree far more with each other than with the other families, so by default
# each family gets one vote, split between its enabled annotators.
FAMILIES = {
    "ScoreGenes": "marker signatures", "AddModuleScore": "marker signatures",
    "CIA": "marker signatures",
    "CellTypist": "reference models", "scANVI": "reference models",
    "scParadise": "reference models",
    "CyteType": "LLM", "CellTypeAI": "LLM",
}
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

    def relation(self, term, ref):
        """How `term` relates to `ref`: same, finer (a descendant of ref),
        coarser (an ancestor of ref) or conflict."""
        if term == ref:
            return "same"
        if ref in self.ancestors(term):
            return "finer"
        if term in self.ancestors(ref):
            return "coarser"
        return "conflict"


# =========================================================================== #
# Calls of each annotator
# =========================================================================== #
def label_columns(method, obs, leiden_col):
    spec = LABEL_COLUMNS[method]
    if isinstance(spec, str):
        return sorted(c for c in obs.columns if c.startswith(spec))
    return [c.format(leiden=leiden_col) for c in spec if c.format(leiden=leiden_col) in obs]


def voter_name(method, col, n_cols):
    if n_cols == 1:
        return method
    return f"{method}:{col.removeprefix(LABEL_COLUMNS[method]).removesuffix('.pkl')}"


def clean_labels(labels: pd.Series, ranked=False) -> pd.Series:
    """Stripped labels, NaN where the cell abstains. For ranked annotators the
    first (best-scoring) label of a comma-separated list."""
    s = labels.astype(str).str.strip()
    if ranked:
        s = s.str.split(",").str[0].str.strip()
    return s.mask(s.str.lower().isin(ABSTAIN))


def read_calls(annotations: dict, cells: pd.Index, leiden_col, mapper: LabelMapper):
    """Per-cell calls of every voter (one per label column).

    Returns (labels, terms, voters, vocabulary): `labels` and `terms` are
    cells x voters frames (NaN = abstain), `voters` maps each voter to its
    annotator, `vocabulary` lists (voter, label, cl_id, mapped_by)."""
    labels, terms, voters, vocab = {}, {}, {}, []
    for method, path in annotations.items():
        obs = ad.read_h5ad(path, backed="r").obs
        common = cells.intersection(obs.index)
        if len(common) < 0.9 * len(cells):
            raise ValueError(f"{method}: only {len(common)}/{len(cells)} cells match "
                             f"the clustered object ({path})")
        cols = label_columns(method, obs, leiden_col)
        if not cols:
            print(f"[consensus] WARNING: {method} has no label column; skipped.")
            continue
        id_col = ID_COLUMNS.get(method, "").format(leiden=leiden_col)
        for col in cols:
            voter = voter_name(method, col, len(cols))
            lab = clean_labels(obs.loc[common, col], ranked=method in RANKED).reindex(cells)
            given = {}
            if id_col in obs:  # annotator-given CL id of each label
                given = (pd.DataFrame({"label": lab.reindex(common),
                                       "id": obs.loc[common, id_col].astype(str)})
                         .dropna().groupby("label")["id"].first().to_dict())
            table = {}
            for label in lab.dropna().unique():
                term, how = mapper(label, given.get(label))
                table[label] = term
                vocab.append((method, voter, label, term, how))
            labels[voter], terms[voter], voters[voter] = lab, lab.map(table), method
    if not voters:
        raise ValueError("No annotator produced labels.")
    vocabulary = pd.DataFrame(vocab, columns=["annotator", "voter", "label", "cl_id", "mapped_by"])
    return pd.DataFrame(labels), pd.DataFrame(terms), voters, vocabulary


def combine_terms(row_terms, mapper: LabelMapper):
    """One annotator's vote for one cell from its voters' terms: the most
    specific term when they lie on one lineage (CellTypist Low "classical
    monocyte" + High "monocyte"), else an equal split between the most
    specific conflicting terms. Returns ((term, weight), ...)."""
    ts = {t for t in row_terms if isinstance(t, str)}
    leaves = sorted(t for t in ts if not any(t != u and t in mapper.ancestors(u) for u in ts))
    return tuple((t, 1.0 / len(leaves)) for t in leaves)


def annotator_cell_votes(terms: pd.DataFrame, voters: dict, mapper: LabelMapper):
    """Long table (cell, annotator, cl_id, w): each annotator's per-cell vote,
    its voters combined with `combine_terms`; abstaining cells have no row."""
    frames = []
    for method in dict.fromkeys(voters.values()):
        sub = terms[[v for v, m in voters.items() if m == method]]
        keys = pd.Series(list(map(tuple, sub.astype(object).where(sub.notna(), None).values)),
                         index=sub.index)
        resolved = {k: combine_terms(k, mapper) for k in keys.unique()}
        rows = [(cell, method, t, w) for cell, k in keys.items() for t, w in resolved[k]]
        frames.append(pd.DataFrame(rows, columns=["cell", "annotator", "cl_id", "w"]))
    return pd.concat(frames, ignore_index=True)


def cluster_votes(cell_votes: pd.DataFrame, clusters: pd.Series):
    """(annotator, cluster, cl_id, share): share of the cluster's cells voting
    for the term; abstaining cells count in the cluster size only."""
    size = clusters.value_counts()
    v = cell_votes.assign(cluster=clusters.reindex(cell_votes["cell"]).values)
    out = v.groupby(["annotator", "cluster", "cl_id"], as_index=False)["w"].sum()
    out["share"] = out["w"] / out["cluster"].map(size)
    return out.drop(columns="w")


def voter_shares(labels: pd.DataFrame, terms: pd.DataFrame, voters: dict, clusters: pd.Series):
    """(annotator, voter, cluster, label, cl_id, share) for every voter."""
    size = clusters.value_counts()
    frames = []
    for voter, method in voters.items():
        df = pd.DataFrame({"cluster": clusters, "label": labels[voter],
                           "cl_id": terms[voter]}).dropna()
        out = df.groupby(["cluster", "label", "cl_id"], as_index=False).size()
        out["share"] = out["size"] / out["cluster"].map(size)
        frames.append(out.drop(columns="size").assign(annotator=method, voter=voter))
    return pd.concat(frames, ignore_index=True)[
        ["annotator", "voter", "cluster", "label", "cl_id", "share"]]


# =========================================================================== #
# The vote
# =========================================================================== #
def cluster_weights(annots, weights=None, groups=None):
    """Weight of each annotator voting in a cluster: its user weight (default
    1) divided by the number of voting annotators of its group (`groups` maps
    annotator -> family; by default every annotator is its own group)."""
    weights, groups = weights or {}, groups or {}
    size = pd.Series([groups.get(a, a) for a in annots]).value_counts()
    return {a: weights.get(a, 1.0) / size[groups.get(a, a)] for a in annots}


def consensus_for_cluster(votes: pd.DataFrame, mapper: LabelMapper, min_agreement=0.5,
                          weights=None, groups=None, labels: pd.DataFrame | None = None):
    """Consensus of one cluster. `votes` has annotator, cl_id and share rows;
    `weights` and `groups` set each annotator's weight (`cluster_weights`);
    `labels` (voter rows with label, cl_id, share) gives the consensus in the
    annotators' own words."""
    annots = list(dict.fromkeys(votes["annotator"]))
    if not annots:
        return {"status": "none", "n_annotators": 0}
    aw = cluster_weights(annots, weights, groups)
    w = votes["annotator"].map(aw) * votes["share"]
    total = sum(aw.values())
    leaf = w.groupby(votes["cl_id"].values).sum()
    support = defaultdict(float)
    for term, s in leaf.items():
        for anc in mapper.ancestors(term):
            support[anc] += s
    passing = [t for t, s in support.items() if s > min_agreement * total + 1e-12]
    top = leaf.idxmax()
    if passing:
        best = max(passing, key=lambda t: (mapper.specificity(t), support[t], t))
        status = "agreed" if best in leaf.index else "partial"
    else:
        best, status = top, "none"
    coarser = sum(s for t, s in leaf.items() if t != best and t in mapper.ancestors(best))
    named = pd.Series(dtype=float)
    if labels is not None and len(labels):
        named = labels[labels["cl_id"] == best].groupby("label")["share"].sum()
    row = {
        "consensus_label": named.idxmax() if len(named) else mapper.name(best),
        "consensus_cl_id": best if best in mapper.cl else "",
        "consensus_cl_name": mapper.name(best),
        "agreement": support[best] / total,
        # votes that don't contradict the consensus (it or a coarser term)
        "compatibility": (support[best] + coarser) / total,
        "status": status,
        "top_label": mapper.name(top),
        "top_label_agreement": leaf[top] / total,
        "n_annotators": len(annots),
    }
    # share of each annotator's cells that support the consensus term
    for a, g in votes.groupby("annotator", sort=False):
        row[f"support:{a}"] = float(g.loc[[best in mapper.ancestors(t) for t in g["cl_id"]],
                                          "share"].sum())
    return row


def _cluster_order(values):
    return sorted(values, key=lambda x: [(0, int(p)) if p.isdigit() else (1, p)
                                         for p in str(x).split("_")])


def consensus_table(votes, shares, mapper, clusters: pd.Series, min_agreement=0.5,
                    stability=None, weights=None, groups=None):
    rows = []
    size = clusters.value_counts()
    for cl in _cluster_order(clusters.unique()):
        row = {"cluster": cl, "n_cells": int(size[cl])}
        if stability is not None:
            row["stability"] = stability.get(cl)
        s = shares[shares["cluster"] == cl]
        row.update(consensus_for_cluster(votes[votes["cluster"] == cl], mapper,
                                         min_agreement, weights, groups, s))
        # each voter's leading label, e.g. "CD4 T (0.93)"
        for voter, g in s.groupby("voter", sort=False):
            r = g.loc[g["share"].idxmax()]
            row[f"call:{voter}"] = f"{r['label']} ({r['share']:.2f})"
        rows.append(row)
    return pd.DataFrame(rows)


# =========================================================================== #
# Diagnostics: annotator agreement and validation against known labels
# =========================================================================== #
def pairwise_agreement(terms: pd.DataFrame, mapper: LabelMapper):
    """For every pair of voters, over the cells both label: the share with the
    same term, and with compatible terms (same, or one an ancestor of the
    other)."""
    rows = []
    for a, b in itertools.combinations(terms.columns, 2):
        pair = terms[[a, b]].dropna()
        if pair.empty:
            continue
        counts = pair.value_counts()
        same = sum(n for (x, y), n in counts.items() if x == y)
        compat = sum(n for (x, y), n in counts.items()
                     if mapper.relation(x, y) != "conflict")
        rows.append((a, b, len(pair), same / len(pair), compat / len(pair)))
    return pd.DataFrame(rows, columns=["voter_a", "voter_b", "n_cells", "same", "compatible"])


OUTCOMES = ["exact", "finer", "coarser", "wrong", "no call"]


def _score(pred: pd.DataFrame, truth: pd.Series, mapper):
    """pred: (cell, cl_id, w) rows; cells missing from pred have no call.
    Returns the weighted outcome counts per cell (cells x OUTCOMES)."""
    rel = {"same": "exact", "finer": "finer", "coarser": "coarser", "conflict": "wrong"}
    out = pd.DataFrame(0.0, index=truth.index, columns=OUTCOMES)
    p = pred[pred["cell"].isin(truth.index)]
    cache = {}
    for cell, t, w in zip(p["cell"], p["cl_id"], p["w"]):
        key = (t, truth[cell])
        if key not in cache:
            cache[key] = rel[mapper.relation(*key)]
        out.at[cell, cache[key]] += w
    out["no call"] = 1.0 - out[OUTCOMES[:-1]].sum(axis=1)
    return out


def validate(truth_labels: pd.Series, cell_votes, cell_consensus: pd.Series,
             clusters: pd.Series, mapper):
    """Score each annotator, the consensus and the best any cluster-level
    labelling could do (each cluster labelled with its majority known label)
    against known labels. `finer` (more specific than the known label, e.g.
    "classical monocyte" for "monocyte") counts as correct."""
    truth_labels = truth_labels.dropna().astype(str)
    truth_labels = truth_labels[~truth_labels.str.lower().isin(ABSTAIN)]
    tmap = {lab: mapper(lab) for lab in truth_labels.unique()}
    truth = truth_labels.map(lambda x: tmap[x][0])
    methods = {a: g[["cell", "cl_id", "w"]] for a, g in cell_votes.groupby("annotator", sort=False)}
    cons = cell_consensus.reindex(truth.index).dropna()
    methods["Consensus"] = pd.DataFrame({"cell": cons.index, "cl_id": cons.values, "w": 1.0})
    best = (pd.DataFrame({"c": clusters.reindex(truth.index), "t": truth})
            .groupby("c")["t"].agg(lambda s: s.value_counts().index[0]))
    ceiling = clusters.reindex(truth.index).map(best)
    methods["Cluster ceiling"] = pd.DataFrame({"cell": ceiling.index, "cl_id": ceiling.values,
                                               "w": 1.0})
    summary, by_label = [], []
    for name, pred in methods.items():
        scored = _score(pred, truth, mapper)
        row = scored.mean()
        summary.append({"method": name, "n_cells": len(truth), **row.to_dict(),
                        "correct": row["exact"] + row["finer"],
                        "correct_or_coarser": row["exact"] + row["finer"] + row["coarser"]})
        per = scored.groupby(truth_labels).mean()
        per["correct"] = per["exact"] + per["finer"]
        for lab, r in per.iterrows():
            by_label.append({"method": name, "known_label": lab,
                             "known_cl_id": tmap[lab][0], "n_cells": int((truth_labels == lab).sum()),
                             **r.to_dict()})
    return pd.DataFrame(summary), pd.DataFrame(by_label)


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
def _name_value(items, what, cast=str):
    out = {}
    for item in items:
        name, sep, value = item.partition("=")
        if not sep or name not in LABEL_COLUMNS:
            raise SystemExit(f"--{what} {item!r}: expected NAME=VALUE with NAME in "
                             f"{list(LABEL_COLUMNS)}")
        out[name] = cast(value)
    return out


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
                   help="Share of the (weighted) voting annotators a term needs more than.")
    p.add_argument("--vote-by", choices=["family", "annotator"], default="family",
                   help="One vote per family of annotators (marker signatures, reference "
                        "models, LLM), split between its members (default), or one per "
                        "annotator.")
    p.add_argument("--weight", action="append", default=[], metavar="NAME=W",
                   help="Multiplies an annotator's vote (default 1).")
    p.add_argument("--truth-col", default="",
                   help="Column of the Step 1 object with known labels: score the "
                        "annotators and the consensus against it.")
    p.add_argument("-o", "--output", required=True, help="Output directory")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    annotations = _name_value(args.annotation, "annotation")
    weights = _name_value(args.weight, "weight", float)
    if any(w < 0 for w in weights.values()):
        raise SystemExit("--weight: weights must be >= 0")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    leiden_col = extract_best_res(args.verdict)
    adata = ad.read_h5ad(args.input)
    clusters = adata.obs[leiden_col].astype(str)
    if args.truth_col and args.truth_col not in adata.obs:
        print(f"[consensus] WARNING: no column {args.truth_col!r} in {args.input}; "
              "no validation against known labels.")
        args.truth_col = ""

    mapper = LabelMapper(CellOntology(args.ontology), args.label_map)
    labels, terms, voters, vocabulary = read_calls(annotations, clusters.index, leiden_col, mapper)
    cell_votes = annotator_cell_votes(terms, voters, mapper)
    votes = cluster_votes(cell_votes, clusters)
    shares = voter_shares(labels, terms, voters, clusters)
    stability = read_stability(Path(args.input).parent, leiden_col)
    annots = list(dict.fromkeys(voters.values()))
    groups = {a: FAMILIES[a] for a in annots} if args.vote_by == "family" else None
    (out / "consensus_settings.json").write_text(json.dumps({
        "resolution": leiden_col, "ontology": str(args.ontology),
        "ontology_version": mapper.cl.version, "n_terms": len(mapper.cl.names),
        "label_maps": [str(p) for p in args.label_map], "min_agreement": args.min_agreement,
        "vote_by": args.vote_by, "annotators": {
            a: {"family": FAMILIES[a], "weight": weights.get(a, 1.0),
                "voters": [v for v, m in voters.items() if m == a], "file": annotations[a]}
            for a in annots},
        "truth_col": args.truth_col or None}, indent=2))
    table = consensus_table(votes, shares, mapper, clusters, args.min_agreement,
                            stability, weights, groups)
    table.insert(1, "resolution", leiden_col)

    # every voter label, its CL term and how it relates to the cluster consensus
    cons = table.set_index("cluster")
    cons_term = cons["consensus_cl_id"].where(cons["consensus_cl_id"] != "",
                                              cons["consensus_cl_name"].radd("free:"))
    shares["cl_name"] = shares["cl_id"].map(mapper.name)
    shares["relation"] = [mapper.relation(t, cons_term[c]) if cons.at[c, "status"] != "none"
                          else "no consensus" for t, c in zip(shares["cl_id"], shares["cluster"])]
    how = vocabulary.set_index(["voter", "label"])["mapped_by"]
    shares["mapped_by"] = [how[(v, lab)] for v, lab in zip(shares["voter"], shares["label"])]
    shares.to_csv(out / "annotation_long.csv", index=False)
    table.to_csv(out / "consensus_by_cluster.csv", index=False)

    n_cells = labels.notna().sum().rename("n_cells")
    counts = (pd.DataFrame({"voter": np.repeat(labels.columns, len(labels)),
                            "label": labels.T.values.ravel()})
              .dropna().value_counts().rename("n_cells").reset_index())
    vocab = vocabulary.merge(counts, on=["voter", "label"], how="left")
    vocab.insert(4, "cl_name", vocab["cl_id"].map(mapper.name))
    (vocab.sort_values(["mapped_by", "voter", "n_cells"], ascending=[False, True, False])
          .to_csv(out / "label_vocabulary.csv", index=False))
    pairwise_agreement(terms, mapper).to_csv(out / "annotator_agreement.csv", index=False)

    status = clusters.map(cons["status"])
    cell_cons = clusters.map(cons_term).where(status != "none")
    for col, key in [("consensus_label", "consensus_label"),
                     ("consensus_cl_id", "consensus_cl_id"),
                     ("consensus_cl_name", "consensus_cl_name"),
                     ("agreement", "consensus_agreement"),
                     ("status", "consensus_status")]:
        adata.obs[key] = clusters.map(cons[col]).values
    # per-cell calls next to the consensus, for plots and the report
    calls = labels.copy()
    calls.insert(0, "cluster", clusters)
    calls["consensus_label"] = adata.obs["consensus_label"]
    calls.to_csv(out / "cell_calls.csv.gz", index_label="cell")
    if args.truth_col:
        summary, by_label = validate(adata.obs[args.truth_col], cell_votes, cell_cons,
                                     clusters, mapper)
        summary.to_csv(out / "validation.csv", index=False)
        by_label.to_csv(out / "validation_by_label.csv", index=False)
        print("[consensus] validation against", args.truth_col)
        print(summary.set_index("method")[["correct", "coarser", "wrong", "no call"]]
              .round(3).to_string())
    adata.write_h5ad(out / "consensus_annotated.h5ad", compression="gzip")

    n_unmapped = vocab.loc[vocab["mapped_by"] == "unmapped", "label"].nunique()
    print(f"[consensus] {leiden_col}: {len(table)} clusters, "
          f"{len(set(voters.values()))} annotators ({len(voters)} voters); "
          f"status {table['status'].value_counts().to_dict()}; "
          f"{n_unmapped} unmapped label(s), see label_vocabulary.csv; "
          f"{int(n_cells.sum())} calls")


if __name__ == "__main__":
    main()
