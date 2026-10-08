#!/usr/bin/env python
"""Behaviour checks for cluster.py. Run through test/run.sh.

One small project runs merge, split, cluster, collect once; every check reads that run.
Junctions start with TGT (Cys) and end with TGG (Trp), like MiXCR's nSeqCDR3.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

CLUSTER = Path(__file__).resolve().parent.parent / "src" / "cluster.py"

# 48nt, so HILARy's cdr3 is 42nt and the 0.2 threshold cuts at 8 mismatches.
BASE = "TGTGCGAGAGGGTTTGACTACTGGAGCGGTATGGACGTCTACTACTGG"
OTHER = "TGTGCAAGAGATCGCAGCAGCTGGTACGGCTTTGACTACTGG"

failures = 0
checks = 0


def ok(what: str, cond: bool) -> None:
    global failures, checks
    checks += 1
    print(f"  {'ok  ' if cond else 'FAIL'} {what}")
    failures += 0 if cond else 1


def mutate(seq: str, *positions: int) -> str:
    out = list(seq)
    for position in positions:
        out[position] = "A" if out[position] != "A" else "C"
    return "".join(out)


def write(path: Path, rows: list[dict]) -> Path:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def stage(*argv: str) -> str:
    done = subprocess.run([sys.executable, str(CLUSTER), *map(str, argv)],
                          capture_output=True, text=True)
    if done.returncode:
        raise SystemExit(f"{argv[0]} failed:\n{done.stdout}\n{done.stderr}")
    return done.stdout


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def clono(key, v, j, junction, **extra):
    return {"sequence_id": key, "v_call": v, "j_call": j, "junction": junction, **extra}


def main(tmp: Path) -> None:
    # Dataset 0: s1 (donor A), s2 (donor B). Keys ending in "g" or "h" test HILARy's id strip.
    ds0 = write(tmp / "ds0.tsv", [
        clono("xigh", "IGHV1-2*02", "IGHJ4*02", BASE),
        clono("g-ig", "IGHV1-2*01", "IGHJ4*02", mutate(BASE, 10, 20)),
        clono("paired", "IGHV1-2*02", "IGHJ4*02", mutate(BASE, 30),
              v_call_light="IGKV1-5*01", j_call_light="IGKJ1*01", junction_light="TGTCAGCAGTGG"),
        clono("vdiff", "IGHV3-23*01", "IGHJ4*02", BASE),
        clono("jdiff", "IGHV1-2*02", "IGHJ6*02", BASE),
        clono("far", "IGHV1-2*02", "IGHJ4*02", mutate(BASE, *range(5, 35, 3))),
        clono("both", "IGHV3S61", "IGHJ4", OTHER),
        clono("zero", "IGHV1-2*02", "IGHJ4*02", BASE),
    ])
    ab0 = write(tmp / "ab0.tsv", [
        {"sample_id": s, "sequence_id": k, "abundance": n} for s, k, n in [
            ("s1", "xigh", 50), ("s1", "g-ig", 30), ("s1", "paired", 10), ("s1", "vdiff", 5),
            ("s1", "jdiff", 5), ("s1", "far", 5), ("s1", "both", 20), ("s1", "zero", 0),
            ("s2", "both", 7)]])
    # Dataset 1: t1 (donor B). Gene names without an allele, which HILARy's split cannot take.
    ds1 = write(tmp / "ds1.tsv", [
        clono("c1", "IGHV3S61", "IGHJ4", mutate(OTHER, 12)),
        clono("c2", "IGHV3S61", "IGHJ4", OTHER + "GGG"),
    ])
    ab1 = write(tmp / "ab1.tsv", [{"sample_id": "t1", "sequence_id": "c1", "abundance": 3},
                                  {"sample_id": "t1", "sequence_id": "c2", "abundance": 1}])
    donors_tsv = write(tmp / "donors.tsv", [{"sample_id": "0_s1", "donor": "A"},
                                            {"sample_id": "0_s2", "donor": "B"},
                                            {"sample_id": "1_t1", "donor": "B"}])
    donors = ["A", "B", "C"]
    donor_args = [arg for d in donors for arg in ("--donor", d)]

    merged, abundance = tmp / "merged.tsv", tmp / "abundance.tsv"
    stage("merge", "--dataset", "0", ds0, ab0, "mixcr", "--dataset", "1", ds1, ab1, "imported",
          "--out-clonotypes", merged, "--out-abundance", abundance,
          "--out-annotations", tmp / "annotations.tsv", "--out-known", tmp / "known.tsv",
          "--out-gene-usage", tmp / "gene-usage.json")
    stage("split", "--clonotypes", merged, "--abundance", abundance, "--donors", donors_tsv,
          "--out-dir", tmp / "split", *donor_args)
    for d in ("lineages", "nodes", "node-links"):
        (tmp / d).mkdir()
    per_donor = {}
    for index, donor in enumerate(donors):
        clones = tmp / f"clones-{index}.tsv"
        stage("cluster", "--clonotypes", tmp / "split" / f"donor-{index}.tsv",
              "--out-clones", clones, "--threshold", "0.2", "--threads", "1",
              "--clone-prefix", f"{donor}/")
        # In bulk mode the tree step hands the heavy clone through as the lineage.
        assigned = read(clones).rename(columns={"clone_id": "lineage_id"})
        assigned["link"] = "1"
        assigned["group_id"] = assigned["sequence_id"]
        assigned.to_csv(tmp / "lineages" / f"donor-{index}.tsv", sep="\t", index=False)
        per_donor[donor] = dict(zip(assigned["sequence_id"], assigned["lineage_id"]))
    out = tmp / "out"
    collect_out = stage("collect", "--annotations", tmp / "annotations.tsv",
          "--lineages-dir", tmp / "lineages", "--nodes-dir", tmp / "nodes",
          "--node-links-dir", tmp / "node-links", "--abundance", abundance,
          "--donors", donors_tsv, "--clonotypes", merged, *donor_args,
          "--dataset", "0", "--dataset", "1", "--per-dataset-dir", out,
          "--out-nodes", tmp / "nodes.tsv", "--out-lineage-stats", tmp / "lineage-stats.tsv",
          "--out-donor-stats", tmp / "donor-stats.json")

    print("== split ==")
    split = [set(read(tmp / "split" / f"donor-{i}.tsv")["sequence_id"]) for i in range(3)]
    ok("a clonotype seen in two donors goes to both", "0_both" in split[0] and "0_both" in split[1])
    ok("one seen in no sample with a donor goes to none", not any("0_zero" in s for s in split))
    ok("a donor with no samples gets an empty table", not split[2])

    print("== cluster ==")
    a, b = per_donor["A"], per_donor["B"]
    ok("every clonotype a donor holds gets a lineage, keys intact despite HILARy's id mangling",
       set(a) == split[0] and set(b) == split[1])
    ok("near junctions share a lineage, whatever the allele or the light chain",
       a["0_xigh"] == a["0_g-ig"] == a["0_paired"])
    ok("another V gene, another J gene or a distant junction splits the lineage",
       len({a["0_xigh"], a["0_vdiff"], a["0_jdiff"], a["0_far"]}) == 4)
    ok("gene names without an allele cluster, across datasets",
       b["0_both"] == b["1_c1"] != b["1_c2"])
    ok("lineage ids carry the donor, so donors cannot collide",
       all(v.startswith("A/") for v in a.values()) and all(v.startswith("B/") for v in b.values()))
    ok("an empty donor clusters to nothing rather than failing", not per_donor["C"])

    print("== collect ==")
    merged_columns = read(merged).columns
    ok("merge keeps dataset names and known antibodies out of the table the early steps read",
       "data_source" not in merged_columns and "is_known" not in merged_columns
       and set(read(tmp / "annotations.tsv")["data_source"]) == {"mixcr", "imported"})
    usage = json.loads((tmp / "gene-usage.json").read_text())
    merged_rows = read(merged)
    first_gene = lambda c: c.split(",")[0].split("*")[0]
    expected = {ds: {} for ds in set(merged_rows["dataset"])}
    for ds, call in zip(merged_rows["dataset"], merged_rows["v_call"]):
        if call:
            expected[ds][first_gene(call)] = expected[ds].get(first_gene(call), 0) + 1
    ok("merge counts each dataset's clonotypes per V gene, alleles dropped",
       {ds: u.get("v", {}) for ds, u in usage.items()} == expected)
    ok("the tree step's known antibodies file carries no dataset names, and no rows without known antibodies",
       list(read(tmp / "known.tsv").columns) == ["sequence_id", "is_known"]
       and read(tmp / "known.tsv").empty)
    ok("collect joins them back: lineages name their data source",
       set(read(tmp / "lineage-stats.tsv")["data_source"]) >= {"mixcr"})
    # The UI reads the last "[==PROGRESS==]" line: a rising percentage, ending at 100%.
    percents = [float(line.rsplit(": ", 1)[1].rstrip("%"))
                for line in collect_out.splitlines() if line.startswith("[==PROGRESS==]")]
    ok("collect reports rising progress, ending at 100%",
       len(percents) > 2 and percents == sorted(percents) and percents[-1] == 100.0)
    members = read(out / "lineages-0.tsv")
    ok("per-dataset tables drop the dataset prefix",
       set(members["sequence_id"]) == {"xigh", "g-ig", "paired", "vdiff", "jdiff", "far", "both"})
    ok("the clonotype in two donors holds a lineage in each",
       set(members.loc[members["sequence_id"] == "both", "lineage_id"]) == {a["0_both"], b["0_both"]})
    stats = read(tmp / "lineage-stats.tsv").set_index("lineage_id")
    ok("lineage stats name the donor and show the id without it",
       stats.loc[a["0_xigh"], "donor"] == "A"
       and stats.loc[a["0_xigh"], "lineage_label"] == a["0_xigh"].removeprefix("A/"))
    ok("cluster size counts the lineage's distinct sequences", stats.loc[a["0_xigh"], "cluster_size"] == "3")
    donor_stats = {d["donor"]: d for d in json.loads((tmp / "donor-stats.json").read_text())}
    ok("donor stats count each donor's clonotypes and lineages",
       [donor_stats[d]["clonotype_count"] for d in donors] == [7, 3, 0]
       and donor_stats["A"]["lineage_count"] == len(set(a.values())))

    # Metadata but nothing placed: the count columns must still be written.
    loader = importlib.util.spec_from_file_location("cluster_under_test", CLUSTER)
    cluster = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(cluster)
    meta = write(tmp / "meta-unplaced.tsv", [{"sample_id": "s1", "timepoint": "d7"}])
    present = pd.DataFrame({"sample_id": ["0_s1"], "sequence_id": ["0_c1"],
                            "lineage_id": ["L1"], "abundance": [0]})
    links = pd.DataFrame({"lineage_id": ["L1"], "sequence_id": ["0_c1"], "node_id": ["1"]})
    got = cluster._node_metadata(argparse.Namespace(sample_metadata=[("0", meta)]), present, links)
    ok("node metadata with nothing placed still carries the count columns",
       list(got.columns) == ["lineage_id", "node_id", "timepoint", "timepoint__count"])

    # Adaptive mode on a donor with no V, J and CDR3-length group of two: HILARy raises a
    # KeyError there, so each clonotype becomes its own lineage and the method says so.
    lone = write(tmp / "lone.tsv", [
        clono("0_a", "IGHV1-2*02", "IGHJ4*02", BASE),
        clono("0_b", "IGHV3-23*01", "IGHJ4*02", OTHER),
        clono("0_c", "IGHV4-34*01", "IGHJ6*02", BASE + "GGG")])
    stage("cluster", "--clonotypes", lone, "--out-clones", tmp / "lone-clones.tsv",
          "--out-method", tmp / "lone-method.json", "--mode", "adaptive", "--threads", "1",
          "--clone-prefix", "D/")
    lone_clones = read(tmp / "lone-clones.tsv")
    ok("adaptive with no group of two: every clonotype its own lineage, under the donor prefix",
       sorted(lone_clones["clone_id"]) == ["D/1", "D/2", "D/3"]
       and json.loads((tmp / "lone-method.json").read_text())["method"] == "singletons")

    # A donor column with no lineage at all: collect must still report the donor, empty.
    empty = tmp / "empty-collect"
    for d in ("lineages", "nodes", "node-links", "out"):
        (empty / d).mkdir(parents=True)
    stage("collect", "--lineages-dir", empty / "lineages", "--nodes-dir", empty / "nodes",
          "--node-links-dir", empty / "node-links",
          "--abundance", write(empty / "abundance.tsv", [{"sample_id": "0_s1", "sequence_id": "0_a", "abundance": 1}]),
          "--donors", write(empty / "donors.tsv", [{"sample_id": "0_s1", "donor": "A"}]),
          "--clonotypes", write(empty / "clonotypes.tsv", [{"sequence_id": "0_a"}]),
          "--donor", "A", "--dataset", "0", "--per-dataset-dir", empty / "out",
          "--out-nodes", empty / "nodes.tsv", "--out-lineage-stats", empty / "stats.tsv",
          "--out-donor-stats", empty / "donor-stats.json")
    ok("collect with a donor column and no lineage reports the donor with nothing in it",
       json.loads((empty / "donor-stats.json").read_text())
       == [{"donor": "A", "clonotype_count": 0, "lineage_count": 0}])

    # A lineage without a tree still lists its members: one parentless node per distinct
    # sequence, carrying the tree step's reason. Tip counts stay those of real trees.
    member = tmp / "member-collect"
    for d in ("lineages", "nodes", "node-links", "builders", "out"):
        (member / d).mkdir(parents=True)
    write(member / "lineages" / "donor-0.tsv", [
        {"sequence_id": "0_a", "lineage_id": "L1", "link": 1, "group_id": "g1"},
        {"sequence_id": "0_b", "lineage_id": "L1", "link": 1, "group_id": "g1"},
        {"sequence_id": "0_c", "lineage_id": "L1", "link": 1, "group_id": "g2"},
        {"sequence_id": "0_d", "lineage_id": "L2", "link": 1, "group_id": "g3"}])
    node_row = {c: "" for c in ["lineage_id", "node_id", "parent_id", "distance", "is_observed",
                                "label", "heavy_sequence", "light_sequence", "node_depth",
                                "terminal_branch_fraction", "parent_descendant_count"]}
    write(member / "nodes" / "donor-0.tsv", [{**node_row, "lineage_id": "L2", "node_id": "1",
                                              "is_observed": "true", "label": "0_d", "node_depth": "1"}])
    write(member / "node-links" / "donor-0.tsv", [{"lineage_id": "L2", "node_id": "1",
                                                   "sequence_id": "0_d", "link": 1, "is_representative": "true"}])
    why = "Its members' alignments differ in length, so they cannot be lined up."
    write(member / "builders" / "donor-0.tsv", [
        {"lineage_id": "L1", "tree_builder": "", "no_tree_reason": why},
        {"lineage_id": "L2", "tree_builder": "pratchet", "no_tree_reason": ""}])
    stage("collect", "--lineages-dir", member / "lineages", "--nodes-dir", member / "nodes",
          "--node-links-dir", member / "node-links", "--builders-dir", member / "builders",
          "--clonotypes", write(member / "clonotypes.tsv", [{"sequence_id": s} for s in ("0_a", "0_b", "0_c", "0_d")]),
          "--dataset", "0", "--per-dataset-dir", member / "out",
          "--out-nodes", member / "nodes.tsv", "--out-lineage-stats", member / "stats.tsv")
    member_nodes = read(member / "nodes.tsv")
    l1 = member_nodes[member_nodes["lineage_id"] == "L1"]
    member_stats = read(member / "stats.tsv").set_index("lineage_id")
    member_links = read(member / "out" / "node-links-0.tsv")
    ok("a lineage without a tree gets one parentless node per distinct sequence, with the reason",
       len(l1) == 2 and (l1["parent_id"] == "").all() and (l1["no_tree_reason"] == why).all())
    ok("its members link to their sequence's node",
       member_links[member_links["lineage_id"] == "L1"].groupby("node_id")["sequence_id"].nunique().sort_values().tolist() == [1, 2])
    ok("the lineage table carries the reason, and only real trees count tips",
       member_stats.loc["L1", "no_tree_reason"] == why and member_stats.loc["L1", "tip_count"] == "0"
       and member_stats.loc["L2", "tip_count"] == "1" and member_stats.loc["L2", "no_tree_reason"] == "")

    # "Nearest known antibody" shows the known antibody's Clone Id where it has one, else its id; the links
    # keep the id, since it keys the known antibody's axis.
    anc = tmp / "known-collect"
    for d in ("lineages", "nodes", "node-links", "known-distances", "out"):
        (anc / d).mkdir(parents=True)
    write(anc / "lineages" / "donor-0.tsv", [
        {"sequence_id": s, "lineage_id": "L1", "link": 1, "group_id": s}
        for s in ("0_a", "0_b", "1_x", "1_y")])
    write(anc / "nodes" / "donor-0.tsv", [{**node_row, "lineage_id": "L1", "node_id": "1",
                                           "is_observed": "true", "label": "0_a", "node_depth": "1"}])
    write(anc / "node-links" / "donor-0.tsv", [{"lineage_id": "L1", "node_id": "1",
                                                "sequence_id": "0_a", "link": 1, "is_representative": "true"}])
    blank = {"known_nt_heavy": 0, "known_aa_light": "", "known_nt_light": ""}
    write(anc / "known-distances" / "donor-0.tsv", [
        {"sequence_id": "0_a", "known_id": "1_x", "known_aa_heavy": 2, **blank},
        {"sequence_id": "0_b", "known_id": "1_y", "known_aa_heavy": 3, **blank}])
    stage("collect", "--lineages-dir", anc / "lineages", "--nodes-dir", anc / "nodes",
          "--node-links-dir", anc / "node-links", "--known-distances-dir", anc / "known-distances",
          "--clonotypes", write(anc / "clonotypes.tsv", [
              {"sequence_id": "0_a", "clone_label": "7"}, {"sequence_id": "0_b", "clone_label": ""},
              {"sequence_id": "1_x", "clone_label": "42"}, {"sequence_id": "1_y", "clone_label": ""}]),
          "--dataset", "0", "--dataset", "1", "--known", "1", "--per-dataset-dir", anc / "out",
          "--out-nodes", anc / "nodes.tsv", "--out-lineage-stats", anc / "stats.tsv")
    nearest = dict(zip(*read(anc / "out" / "known-distances-0.tsv")[["sequence_id", "known_id"]].T.values))
    linked = sorted(read(anc / "out" / "known-links-0-1.tsv")["known_id"])
    ok("nearest known antibody shows the known antibody's Clone Id, else its id; links keep the id",
       nearest == {"a": "42", "b": "y"} and linked == ["x", "y"])

    # The full method compares rows base by base, so padding must not read as a difference:
    # a member with shorter 5' and 3' coverage pads with its gene's germline, not N.
    v_side, j_side = "ACGTACGTACGTAAAACCCCGGGG", "TTTGGGCCC"
    junction_n = "N" * 12
    def row(sid, cut_v, cut_j):
        g = v_side[cut_v:] + junction_n + j_side[:len(j_side) - cut_j]
        return {"sequence_id": sid, "locus": "IGH", "v_call": "IGHV1-2*02", "j_call": "IGHJ4*02",
                "germline_alignment": g, "masked_sequence_alignment": g,
                "frame_left": str(len(v_side) - cut_v), "frame_right": str(len(j_side) - cut_j)}
    aligned_rows = pd.DataFrame([row("0_full", 0, 0), row("0_short", 8, 4)])
    alt = cluster._alt_alignments(pd.DataFrame({"sequence_id": ["0_full", "0_short"]}), aligned_rows)
    s_full, s_short = alt["alt_sequence_alignment"].tolist()
    ok("full method padding: unmutated members of different coverage are identical, not apart by their padding",
       len(s_full) == len(s_short) and s_full == s_short
       and alt["alt_sequence_alignment"].tolist() == alt["alt_germline_alignment"].tolist())

    # A node's datasets, from merge's data_source, in the pass that counts its clonotypes.
    placed = pd.DataFrame({"sequence_id": ["0_a", "1_b"], "lineage_id": ["L1", "L1"],
                           "abundance": [2.0, 3.0]})
    node_links = pd.DataFrame({"lineage_id": ["L1", "L1", "L1"], "node_id": ["1", "1", "2"],
                               "sequence_id": ["0_a", "1_b", "0_a"]})
    sources = pd.DataFrame({"sequence_id": ["0_a", "1_b"], "data_source": ["mixcr", "imported"]})
    got = cluster._node_abundance(placed, node_links, sources).set_index("node_id")
    ok("a node seen in two datasets names both, sorted; one seen in one names it",
       got.loc["1", "dataset"] == "imported, mixcr" and got.loc["2", "dataset"] == "mixcr")
    ok("and its clonotype count and abundance are unchanged",
       got.loc["1", "clonotype_count"] == 2 and got.loc["1", "abundance"] == 5.0)

    # HILARy's pool pickles the seeded simulation by name; a worker must find the seeded one,
    # or every worker dies unpickling and the pool respawns them forever.
    import pickle
    namespace = {}
    exec(cluster.HILARY_SEEDING, namespace)
    from hilary import inference
    method = object.__new__(inference.HILARy).simulate_xs_ys
    back = pickle.loads(pickle.dumps(method))
    ok("the seeded HILARy simulation survives the pool's pickling",
       getattr(back, "__wrapped__", None) is namespace["simulate"])

    # Run id: row order must not change it, membership must.
    members = pd.DataFrame({"sequence_id": ["a", "b"], "lineage_id": ["L1", "L1"]})
    tree = pd.DataFrame({"lineage_id": ["L1"], "node_id": ["1"], "parent_id": [""], "label": ["a"]})
    moved = members.assign(lineage_id=["L1", "L2"])
    # Observed nodes take their representative's short clone label; others keep theirs.
    nodes = pd.DataFrame({"lineage_id": ["L1", "L1", "L1"], "node_id": ["1", "2", "3"],
                          "label": ["key-a", "key-b", ""]})
    node_links = pd.DataFrame({"lineage_id": ["L1", "L1", "L1"], "node_id": ["1", "1", "2"],
                               "sequence_id": ["0_a", "0_a2", "0_b"],
                               "is_representative": ["true", "false", "true"]})
    clones = pd.DataFrame({"sequence_id": ["0_a", "0_a2", "0_b"],
                           "clone_label": ["C-ABCDE", "C-OTHER", ""]})
    ok("observed nodes are named by their clone label, others keep their label",
       list(cluster._short_labels(nodes, node_links, clones)) == ["C-ABCDE", "key-b", ""])
    ok("run id ignores row order and follows lineage membership",
       cluster._run_id(members, tree) == cluster._run_id(members.iloc[::-1], tree)
       and cluster._run_id(members, tree) != cluster._run_id(moved, tree))


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        main(Path(tmp))
    print(f"\n{checks} checks, {failures} failures")
    sys.exit(1 if failures else 0)
