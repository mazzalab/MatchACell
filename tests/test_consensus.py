"""Step 3 consensus: ontology reading, label mapping and hierarchical voting."""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
import consensus_tool as ct  # noqa: E402
from cell_ontology import CellOntology  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

OBO = """format-version: 1.2

[Term]
id: CL:0000000
name: cell

[Term]
id: CL:0000738
name: leukocyte
is_a: CL:0000000 ! cell

[Term]
id: CL:0000084
name: T cell
synonym: "T lymphocyte" EXACT []
synonym: "thymocyte-like" RELATED []
synonym: "TC2" EXACT []
is_a: CL:0000738 ! leukocyte

[Term]
id: CL:0000624
name: CD4-positive, alpha-beta T cell
is_a: CL:0000084 ! T cell

[Term]
id: CL:0000625
name: CD8-positive, alpha-beta T cell
is_a: CL:0000084 ! T cell

[Term]
id: CL:0000236
name: B cell
is_a: CL:0000738 ! leukocyte

[Term]
id: CL:0000999
name: retired cell
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


@pytest.fixture
def obo(tmp_path):
    path = tmp_path / "cl.obo"
    path.write_text(OBO)
    return path


@pytest.fixture
def mapper(obo, tmp_path):
    label_map = tmp_path / "map.tsv"
    label_map.write_text("label\tcl_id\nCD4 T\tCL:0000624\nCD8 T\tCL:0000625\n")
    return ct.LabelMapper(CellOntology(obo), [label_map])


def test_ontology_reads_is_a_names_and_exact_synonyms(obo):
    cl = CellOntology(obo)
    assert "CL:0000999" not in cl                       # obsolete
    assert cl.ancestors("CL:0000624") == {"CL:0000624", "CL:0000084", "CL:0000738", "CL:0000000"}
    assert cl.specificity("CL:0000624") > cl.specificity("CL:0000084")
    assert cl.lookup("T Lymphocytes") == "CL:0000084"     # synonym, case, plural
    assert cl.lookup("B") == "CL:0000236"                 # "<label> cell"
    assert cl.lookup("thymocyte-like") is None            # RELATED synonyms don't count
    assert cl.lookup("TC2") is None                       # abbreviations: names only


def test_mapper_order_and_unknown_ids(mapper, obo, tmp_path):
    assert mapper("cd4 t") == ("CL:0000624", "label map")
    assert mapper("B cells") == ("CL:0000236", "CL name")
    assert mapper("whatever", cl_id="CL:0000625") == ("CL:0000625", "annotator")
    assert mapper("Mystery blob") == ("free:mystery blob", "unmapped")
    bad = tmp_path / "bad.tsv"
    bad.write_text("label\tcl_id\nX\tCL:1234567\n")
    with pytest.raises(ValueError, match="CL:1234567"):
        ct.LabelMapper(CellOntology(obo), [bad])


def test_ranked_labels_keep_the_best_and_abstentions_drop():
    labels = pd.Series(["CD4 T, CD8 T", "CD4 T", "Unknown", " B cells "])
    assert ct.clean_labels(labels, ranked=True).tolist()[:2] == ["CD4 T", "CD4 T"]
    assert ct.clean_labels(labels).isna().tolist() == [False, False, True, False]
    assert ct.clean_labels(labels).iloc[3] == "B cells"


def test_models_of_one_annotator_combine_along_the_lineage(mapper):
    # CellTypist Low + High: the finer term when they agree, a split when not
    assert ct.combine_terms(("CL:0000624", "CL:0000084"), mapper) == (("CL:0000624", 1.0),)
    assert ct.combine_terms(("CL:0000624", "CL:0000236"), mapper) == (
        ("CL:0000236", 0.5), ("CL:0000624", 0.5))
    assert ct.combine_terms((None, None), mapper) == ()


def test_cluster_votes_count_abstaining_cells_in_the_cluster_size():
    clusters = pd.Series(["0", "0", "0", "1"], index=list("abcd"))
    cell_votes = pd.DataFrame({"cell": ["a", "b", "b", "d"], "annotator": "x",
                               "cl_id": ["T", "T", "B", "B"], "w": [1.0, 0.5, 0.5, 1.0]})
    v = ct.cluster_votes(cell_votes, clusters).set_index(["cluster", "cl_id"])["share"]
    assert v[("0", "T")] == pytest.approx(1.5 / 3)
    assert v[("0", "B")] == pytest.approx(0.5 / 3)
    assert v[("1", "B")] == 1.0


def _votes(rows):
    return pd.DataFrame(rows, columns=["annotator", "cl_id", "label", "share"])


def test_consensus_takes_the_most_specific_majority_term(mapper):
    v = _votes([("a", "CL:0000624", "CD4 T", 1.0), ("b", "CL:0000624", "CD4 T", 0.9),
                ("c", "CL:0000084", "T cells", 1.0)])
    r = ct.consensus_for_cluster(v, mapper, labels=v)
    assert (r["consensus_cl_id"], r["status"], r["consensus_label"]) == ("CL:0000624", "agreed", "CD4 T")
    assert r["agreement"] == pytest.approx(1.9 / 3)
    assert r["compatibility"] == pytest.approx(1.0 - 0.1 / 3)  # the coarse T cell vote fits
    assert r["support:a"] == 1.0 and r["support:c"] == 0.0


def test_consensus_falls_back_to_the_common_ancestor(mapper):
    v = _votes([("a", "CL:0000624", "CD4 T", 1.0), ("b", "CL:0000625", "CD8 T", 1.0),
                ("c", "CL:0000236", "B cells", 1.0)])
    r = ct.consensus_for_cluster(v, mapper)
    assert r["consensus_cl_name"] == "T cell" and r["status"] == "partial"
    assert r["agreement"] == pytest.approx(2 / 3)


def test_family_votes_split_between_the_members_that_voted(mapper):
    # two marker methods agree on CD4, the one reference model says CD8:
    # per annotator CD4 wins 2:1, per family it is 1:1 and only T cell passes
    v = _votes([("ScoreGenes", "CL:0000624", "CD4 T", 1.0),
                ("AddModuleScore", "CL:0000624", "CD4 T", 1.0),
                ("CellTypist", "CL:0000625", "CD8 T", 1.0)])
    assert ct.consensus_for_cluster(v, mapper)["consensus_cl_name"].startswith("CD4")
    r = ct.consensus_for_cluster(v, mapper, groups=ct.FAMILIES)
    assert (r["consensus_cl_name"], r["status"]) == ("T cell", "partial")
    assert ct.cluster_weights(["ScoreGenes", "CellTypist"], {"CellTypist": 2.0}, ct.FAMILIES) == {
        "ScoreGenes": 1.0, "CellTypist": 2.0}


def test_consensus_without_a_majority(mapper):
    v = _votes([("a", "free:blob", "blob", 1.0), ("b", "free:goo", "goo", 1.0)])
    r = ct.consensus_for_cluster(v, mapper)
    assert r["status"] == "none" and r["consensus_cl_id"] == ""


def test_relation_and_pairwise_agreement(mapper):
    assert mapper.relation("CL:0000624", "CL:0000084") == "finer"
    assert mapper.relation("CL:0000084", "CL:0000624") == "coarser"
    assert mapper.relation("CL:0000624", "CL:0000625") == "conflict"
    terms = pd.DataFrame({"x": ["CL:0000624", "CL:0000624", "CL:0000236", None],
                          "y": ["CL:0000624", "CL:0000084", "CL:0000624", "CL:0000236"]})
    r = ct.pairwise_agreement(terms, mapper).iloc[0]
    assert (r["n_cells"], r["same"], r["compatible"]) == (3, pytest.approx(1 / 3), pytest.approx(2 / 3))


def _write(tmp_path, name, obs):
    a = ad.AnnData(np.zeros((len(obs), 2), dtype="float32"), obs=obs)
    path = tmp_path / f"{name}.h5ad"
    a.write_h5ad(path)
    return path


def test_main_end_to_end(obo, tmp_path):
    cells = [f"c{i}" for i in range(6)]
    clusters = pd.Series(["0", "0", "0", "1", "1", "1"], index=cells)
    clustered = _write(tmp_path, "clustered_multi_resolution", pd.DataFrame(
        {"leiden_0.5": clusters, "known": ["CD4 T"] * 3 + ["B cells"] * 2 + ["CD8 T"]}))
    verdict = tmp_path / "MatchA_Verdict.txt"
    verdict.write_text("Recommended resolution : leiden_0.5\n")
    sg = _write(tmp_path, "sg", pd.DataFrame(
        {"cell_type_pred": ["CD4 T, CD8 T"] * 3 + ["B cells"] * 3}, index=cells))
    ctyp = _write(tmp_path, "ct", pd.DataFrame({
        "celltypist_label_majority_voting_Low.pkl": ["CD4 T"] * 3 + ["B cells"] * 3,
        "celltypist_label_majority_voting_High.pkl": ["T cells"] * 3 + ["B cells"] * 3,
    }, index=cells))
    cia = _write(tmp_path, "cia", pd.DataFrame(
        {"CIA prediction default": ["CD4 T", "CD4 T", "Unknown", "B cells", "B cells", "Odd"]},
        index=cells))
    label_map = tmp_path / "map.tsv"
    label_map.write_text("label\tcl_id\nCD4 T\tCL:0000624\nCD8 T\tCL:0000625\n")
    out = tmp_path / "out"

    ct.main(["-i", str(clustered), "-v", str(verdict), "--ontology", str(obo),
             "--label-map", str(label_map), "-o", str(out), "--truth-col", "known",
             "-a", f"ScoreGenes={sg}", "-a", f"CellTypist={ctyp}", "-a", f"CIA={cia}"])

    table = pd.read_csv(out / "consensus_by_cluster.csv", dtype={"cluster": str}).set_index("cluster")
    assert table.loc["0", "consensus_label"] == "CD4 T"
    assert table.loc["1", "consensus_cl_name"] == "B cell" and table.loc["1", "status"] == "agreed"
    assert {"call:ScoreGenes", "call:CellTypist:Low", "call:CIA"} <= set(table.columns)
    vocab = pd.read_csv(out / "label_vocabulary.csv")
    assert vocab.loc[vocab["mapped_by"] == "unmapped", "label"].tolist() == ["Odd"]
    assert set(pd.read_csv(out / "annotation_long.csv")["relation"]) <= {
        "same", "finer", "coarser", "conflict"}
    val = pd.read_csv(out / "validation.csv").set_index("method")
    assert val.loc["Consensus", "correct"] == pytest.approx(5 / 6)   # the CD8 cell is wrong
    assert val.loc["CIA", "no call"] == pytest.approx(1 / 6)
    assert val.loc["Cluster ceiling", "correct"] == pytest.approx(5 / 6)
    calls = pd.read_csv(out / "cell_calls.csv.gz", index_col="cell")
    assert calls.loc["c0", "ScoreGenes"] == "CD4 T" and calls.loc["c0", "consensus_label"] == "CD4 T"
    assert (out / "annotator_agreement.csv").exists()
    assert json.loads((out / "consensus_settings.json").read_text())["vote_by"] == "family"
    obs = ad.read_h5ad(out / "consensus_annotated.h5ad").obs
    assert obs.loc["c4", "consensus_label"] == "B cells"
    assert {"consensus_cl_id", "consensus_agreement", "consensus_status"} <= set(obs.columns)

    # the Steps 1-3 report reads the same folders
    import report_tool
    a = ad.read_h5ad(clustered)
    a.obsm["X_umap"] = np.arange(12, dtype=float).reshape(6, 2)
    a.write_h5ad(clustered)
    page = tmp_path / "report.html"
    report_tool.main(["--sample", "s1", "--step1", str(tmp_path), "--consensus", str(out),
                      "--settings", '{"type": "raw"}', "-o", str(page)])
    text = page.read_text()
    for part in ("Step 1", "Step 2", "Step 3", "Calls per cluster", "Against known labels",
                 "MATCHA_EMBED", "CL:0000624"):
        assert part in text, part


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="Snakemake is not installed")
def test_consensus_target_schedules_enabled_annotators(tmp_path):
    import yaml

    config = yaml.safe_load((ROOT / "config.yaml").read_text())
    config["output_dir"] = str(tmp_path / "out dir")
    sample = tmp_path / "input.h5ad"
    sample.touch()
    config["samples"] = {"s1": str(sample)}
    config["matchacell_annotation"]["annot_file"] = str(sample)
    config["matchacell_annotation"]["celltypeai"]["tissue"] = ""
    config["matchacell_consensus"]["exclude"] = ["CIA"]
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    env = {**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")}
    env.pop("CYTETYPE_API_TOKEN", None)
    result = subprocess.run(
        ["snakemake", "--snakefile", str(ROOT / "Snakefile"), "--configfile", str(path),
         "--cores", "1", "--dry-run", "--printshellcmds", "consensus"],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r"^matchacell_consensus\s+1", result.stdout, re.MULTILINE)
    assert re.search(r"^fetch_cell_ontology\s+1", result.stdout, re.MULTILINE)
    cmd = result.stdout[result.stdout.index("consensus_tool.py"):]
    assert "ScoreGenes=" in cmd and "CellTypist=" in cmd
    assert "CIA=" not in cmd and "CellTypeAI=" not in cmd


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="Snakemake is not installed")
@pytest.mark.parametrize("method", ["scDblFinder", ""])
def test_doublet_rule_runs_only_when_enabled(tmp_path, method):
    import yaml

    config = yaml.safe_load((ROOT / "config.yaml").read_text())
    config["output_dir"] = str(tmp_path / "out")
    sample = tmp_path / "input.h5ad"
    sample.touch()
    config["samples"] = {"s1": str(sample)}
    config["matchacell_cluster_stability"]["doublets"].update(
        {"method": method, "samples_col": "lane", "action": "flag"})
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    env = {**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")}
    result = subprocess.run(
        ["snakemake", "--snakefile", str(ROOT / "Snakefile"), "--configfile", str(path),
         "--cores", "1", "--dry-run", "--printshellcmds", "cluster_stability"],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    scheduled = re.search(r"^doublets\s+1", result.stdout, re.MULTILINE) is not None
    assert scheduled == bool(method)
    step1 = result.stdout[result.stdout.index("matchacell_cluster_stability.py"):]
    if method:
        assert "scdblfinder_tool.R" in result.stdout and "--samples-col lane" in result.stdout
        assert "--doublets" in step1 and "--doublet-action flag" in step1
    else:
        assert "--doublets" not in step1
