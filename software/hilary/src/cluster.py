#!/usr/bin/env python
"""Heavy-chain clone assignment with HILARy, and the per-lineage summaries.

Stages, in workflow order:

  merge    per-dataset clonotype and abundance tables -> one of each.
           Ids get the dataset's position as a prefix, so they cannot collide.
  split    clonotypes + abundance + donors -> one clonotype table per donor.
           A lineage never spans donors; abundance says whose a clonotype is.
  cluster  clonotypes (sequence_id, v_call, j_call, junction)
           -> --out-clones (sequence_id, clone_id).
           fixed: crude method at --threshold. adaptive: full method at
           --precision/--sensitivity with --aligned, or cdr3 method if some
           clonotype has no alignment. --out-method records which ran.
  collect  per-donor tree outputs + abundance -> merged tables:
             --out-nodes            every donor's nodes
             --out-lineage-stats    lineage_id, cluster_size, tip_count, tree_builder, no_tree_reason,
                                    genes, CDR3, abundance, sample count, donor, source
             --out-node-properties  lineage_id, node_id, clonotype content for the tooltip
             --out-node-metadata    lineage_id, node_id, one column per metadata
                                    column, values joined with ", "
             --out-run-id           JSON content id of the lineages and trees
             --out-donor-stats      JSON clonotype and lineage counts per donor
           and per dataset under --per-dataset-dir, ids unprefixed:
             lineages-<i>.tsv            sequence_id, lineage_id, link
             node-links-<i>.tsv          node-to-clonotype linker
             germline-mutations-<i>.tsv  sequence_id, germline_mutation_count
             member-counts-<i>.tsv       lineage_id, member_count
             known-distances-<i>.tsv    sequence_id, known_id, per-chain aa and nt
                                         mutations to that known antibody

Clustering uses heavy chains only, so bulk and single-cell clonotypes can share
a lineage. Light chains are applied later by Dowser in `trees.R`.
The adaptive methods are calibrated on human data only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

# HILARy's own name for the clustering it writes back.
HILARY_CLONE_COLUMN = "clone_id"

# What a progress line starts with; the block shows the last one in its bar.
PROGRESS_PREFIX = "[==PROGRESS==]"


STARTED = time.monotonic()


def clock() -> str:
    """Elapsed wall time since this run started, as H:MM:SS."""
    s = round(time.monotonic() - STARTED)
    return f"{s // 3600}:{s // 60 % 60:02d}:{s % 60:02d}"


# Time goes first: the UI reads a percentage off the end.
def progress(text: str) -> None:
    print(f"{PROGRESS_PREFIX} [{clock()}] {text}", flush=True)

# `junction`, not `cdr3`: MiXCR's CDR3 includes Cys and Trp/Phe (AIRR's junction),
# and HILARy derives cdr3 itself as junction[3:-3].
REQUIRED_COLUMNS = ("sequence_id", "v_call", "j_call", "junction")


def _write_hilary_input(frame: pd.DataFrame, path: Path) -> None:
    """Write HILARy's input, with placeholder alignments if none are real.

    HILARy always computes a mutation count and raises without alignment columns;
    only the full method reads it, and that runs only with real alignments.
    """
    prepared = frame.copy()

    # HILARy's own allele split raises when no call has a "*" allele suffix,
    # which is normal for our gene-level input.
    for call, gene in (("v_call", "v_gene"), ("j_call", "j_gene")):
        if gene not in prepared.columns:
            prepared[gene] = prepared[call].str.split("*", n=1).str[0]

    # HILARy prefers `alt_` columns, so adding them always would hide real alignments.
    real = ("v_sequence_alignment", "j_sequence_alignment",
            "v_germline_alignment", "j_germline_alignment")
    has_alt = "alt_sequence_alignment" in prepared.columns
    has_real = all(column in prepared.columns for column in real)

    if not has_alt and not has_real:
        prepared["alt_sequence_alignment"] = prepared["junction"]
        prepared["alt_germline_alignment"] = prepared["junction"]

    prepared.to_csv(path, sep="\t", index=False)


# HILARy's command per method, and the file each writes next to the input's name.
HILARY_COMMANDS = {
    "crude": ("crude-method", "inferred_crude_method_"),
    "cdr3": ("cdr3-method", "inferred_cdr3_based_"),
    "full": ("full-method", "inferred_full_method_"),
}


# HILARy has no seed option. Its adaptive methods subsample large classes in the parent and
# simulate null distributions per class in forked workers, which share the parent's state and
# take classes in whatever order they free up. So the parent is seeded, and each class's
# simulation is seeded by its class id.
HILARY_SEED = 0
HILARY_SEEDING = f"""
import functools
import zlib
import numpy as np
from hilary import inference

simulate = inference.HILARy.simulate_xs_ys

# Wrapped, so the pool pickles it by its own name and workers find the seeded one.
@functools.wraps(simulate)
def seeded(self, args):
    np.random.seed(({HILARY_SEED} + zlib.crc32(str(args[-1]).encode())) % 2**32)
    return simulate(self, args)

inference.HILARy.simulate_xs_ys = seeded
np.random.seed({HILARY_SEED})
"""
HILARY_BOOTSTRAP = HILARY_SEEDING + """
from hilary.__main__ import app
app(prog_name="hilary")
"""


def run_hilary(
    input_tsv: Path,
    work_dir: Path,
    method: str,
    threshold: float,
    precision: float,
    sensitivity: float,
    threads: int,
) -> Path:
    """Run one HILARy method; return the path of its output.

    The crude threshold is a fraction of CDR3 length, not a count.
    """
    command_name, output_prefix = HILARY_COMMANDS[method]
    command = [
        sys.executable,
        "-c",
        HILARY_BOOTSTRAP,
        command_name,
        str(input_tsv),
        "--result-folder",
        str(work_dir),
        "--threads",
        str(threads),
        "--override",
    ]
    if method == "crude":
        # Underscore, not hyphen: HILARy's own spelling.
        command += ["--normalized_threshold", str(threshold)]
    else:
        # `--silent` exists on the adaptive commands alone.
        command += ["--precision", str(precision), "--sensitivity", str(sensitivity), "--silent"]
    subprocess.run(command, check=True)

    written = work_dir / f"{output_prefix}{input_tsv.name}"
    if not written.exists():
        raise RuntimeError(f"HILARy reported success but {written} is missing")
    return written


def _alt_alignments(clonotypes: pd.DataFrame, aligned: pd.DataFrame) -> pd.DataFrame | None:
    """Heavy alignments for the full method, padded to one length; None if any are missing.

    HILARy compares rows position by position, so they must line up.
    """
    if aligned.empty or "masked_sequence_alignment" not in aligned.columns:
        return None
    heavy = aligned[aligned["locus"] == "IGH"].drop_duplicates("sequence_id").set_index("sequence_id")
    if not clonotypes["sequence_id"].isin(heavy.index).all():
        return None
    left = heavy["frame_left"].astype(int)
    right = heavy["frame_right"].astype(int)
    germline = heavy["germline_alignment"].astype(str)
    # Rows end at different 5' and 3' points. HILARy compares two rows base by base, and an N
    # against a base counts as a difference, so N padding reads as mutations and keeps related
    # members of different coverage apart. A row is padded with its gene's germline instead,
    # the same in sequence and germline, so the padding reads as unmutated. Each flank is taken
    # from the longest germline of that gene, lined up at the junction; N only where no row of
    # the gene reaches.
    gene = lambda calls: calls.astype(str).str.split(",", n=1).str[0].str.split("*", n=1).str[0]
    v_gene, j_gene = gene(heavy["v_call"]), gene(heavy["j_call"])
    v_flank = pd.Series([g[:l] for g, l in zip(germline, left)], index=heavy.index)
    j_flank = pd.Series([g[len(g) - r:] if r else "" for g, r in zip(germline, right)], index=heavy.index)
    longest = lambda flanks, genes: flanks.groupby(genes).agg(lambda f: max(f, key=len)).to_dict()
    v_ref, j_ref = longest(v_flank, v_gene), longest(j_flank, j_gene)
    max_left, max_right = int(left.max()), int(right.max())

    def five_prime(ref: str, have: int) -> str:
        need = max_left - have
        return ("N" * need + ref[:len(ref) - have])[-need:] if need > 0 else ""

    def three_prime(ref: str, have: int) -> str:
        need = max_right - have
        return (ref[have:] + "N" * need)[:need] if need > 0 else ""

    pad_left = pd.Series([five_prime(v_ref[g], l) for g, l in zip(v_gene, left)], index=heavy.index)
    pad_right = pd.Series([three_prime(j_ref[g], r) for g, r in zip(j_gene, right)], index=heavy.index)
    padded = pd.DataFrame({
        "alt_sequence_alignment": pad_left + heavy["masked_sequence_alignment"] + pad_right,
        "alt_germline_alignment": pad_left + germline + pad_right,
    })
    return padded.reindex(clonotypes["sequence_id"]).set_index(clonotypes.index)


# Tables are read with `keep_default_na=False`: a donor called "NA" is a donor.

# Between a dataset's position and an id. Not ":", which Newick reserves in tip labels.
DATASET_SEP = "_"


def prefixed(index: str, values: pd.Series) -> pd.Series:
    return index + DATASET_SEP + values.astype(str)


def dataset_of(values: pd.Series) -> pd.Series:
    """The dataset position an id was prefixed with."""
    return values.astype(str).str.split(DATASET_SEP, n=1).str[0]


def unprefixed(values: pd.Series) -> pd.Series:
    """The id as the dataset knows it."""
    return values.astype(str).str.split(DATASET_SEP, n=1).str[1]


def merge(args: argparse.Namespace) -> None:
    """Merge every dataset into one clonotype and one abundance table; missing columns read empty."""
    clonotype_parts, abundance_parts = [], []
    for index, clonotypes_path, abundance_path, data_source in args.dataset:
        clonotypes = pd.read_csv(clonotypes_path, sep="\t", dtype=str, keep_default_na=False)
        clonotypes["sequence_id"] = prefixed(index, clonotypes["sequence_id"])
        clonotypes["dataset"] = index
        clonotypes["data_source"] = data_source
        # Known antibodies are known antibodies; the tree step keeps them and measures mates to them.
        clonotypes["is_known"] = "true" if index in args.known else "false"
        clonotype_parts.append(clonotypes)
        abundance = pd.read_csv(abundance_path, sep="\t", dtype=str, keep_default_na=False)
        abundance["sequence_id"] = prefixed(index, abundance["sequence_id"])
        abundance["sample_id"] = prefixed(index, abundance["sample_id"])
        abundance_parts.append(abundance)
        print(f"dataset {index}: {len(clonotypes)} clonotypes, {len(abundance)} abundance rows "
              f"({data_source})", file=sys.stderr)
    merged = pd.concat(clonotype_parts, ignore_index=True).fillna("")
    # Clonotypes per heavy V and J gene and dataset, so the block can tell whether two datasets
    # name genes alike: clustering only compares clonotypes of the same V and J gene.
    if args.out_gene_usage is not None:
        first_gene = lambda calls: calls.astype(str).str.split(",", n=1).str[0].str.split("*", n=1).str[0]
        usage: dict = {}
        for kind in ("v", "j"):
            column = f"{kind}_call"
            if column not in merged.columns:
                continue
            genes = merged.assign(gene=first_gene(merged[column]))
            for (dataset, gene), count in genes[genes["gene"] != ""].groupby(["dataset", "gene"]).size().items():
                usage.setdefault(str(dataset), {}).setdefault(kind, {})[gene] = int(count)
        args.out_gene_usage.write_text(json.dumps(usage, sort_keys=True))
    # Dataset names and known antibodies go apart: only the tree step and collect read them, so
    # renaming a dataset or toggling a known antibody leaves the table alignment, allele inference
    # and clustering read unchanged, and their results are reused.
    annotation_columns = ["sequence_id", "data_source", "is_known"]
    if args.out_annotations is not None:
        merged[annotation_columns].to_csv(args.out_annotations, sep="\t", index=False)
        # The tree step reads known antibodies only, so renaming a dataset leaves its input unchanged.
        if args.out_known is not None:
            merged.loc[merged["is_known"] == "true", ["sequence_id", "is_known"]].to_csv(
                args.out_known, sep="\t", index=False)
        merged = merged.drop(columns=["data_source", "is_known"])
    merged.to_csv(args.out_clonotypes, sep="\t", index=False)
    pd.concat(abundance_parts, ignore_index=True).fillna("").to_csv(
        args.out_abundance, sep="\t", index=False
    )


def split(args: argparse.Namespace) -> None:
    """Write donor-<i>.tsv per donor, by position; a clonotype seen in several donors goes to each."""
    clonotypes = pd.read_csv(args.clonotypes, sep="\t", dtype=str, keep_default_na=False)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    def write(index: int, frame: pd.DataFrame) -> None:
        frame.to_csv(args.out_dir / f"donor-{index}.tsv", sep="\t", index=False)

    # No donor column: everything is one donor.
    if not args.donor:
        write(0, clonotypes)
        print(f"single donor: {len(clonotypes)} clonotypes", file=sys.stderr)
        return

    abundance = pd.read_csv(args.abundance, sep="\t", dtype={"sample_id": str, "sequence_id": str},
                            keep_default_na=False)
    abundance = abundance[pd.to_numeric(abundance["abundance"], errors="coerce").fillna(0) > 0]
    donors = pd.read_csv(args.donors, sep="\t", dtype=str, keep_default_na=False)
    abundance = abundance.merge(donors, on="sample_id", how="inner")

    seen = abundance[["donor", "sequence_id"]].drop_duplicates()
    per_donor = seen.groupby("sequence_id")["donor"].nunique()
    shared = int((per_donor > 1).sum())
    placed = set(seen["sequence_id"])
    unplaced = len(set(clonotypes["sequence_id"]) - placed)

    for index, donor in enumerate(args.donor):
        ids = set(seen.loc[seen["donor"] == donor, "sequence_id"])
        subset = clonotypes[clonotypes["sequence_id"].isin(ids)]
        write(index, subset)
        print(f"donor {donor!r}: {len(subset)} clonotypes", file=sys.stderr)

    if shared:
        print(
            f"{shared} clonotypes were seen in more than one donor and are clustered "
            f"in each of them",
            file=sys.stderr,
        )
    if unplaced:
        print(
            f"{unplaced} clonotypes have no abundance in any sample with a donor and "
            f"are left out",
            file=sys.stderr,
        )


def cluster(args: argparse.Namespace) -> None:
    """Assign every clonotype a heavy-chain clone id."""
    clonotypes = pd.read_csv(args.clonotypes, sep="\t", dtype=str, keep_default_na=False)
    missing = [column for column in REQUIRED_COLUMNS if column not in clonotypes.columns]
    if missing:
        raise SystemExit(f"input is missing required columns: {', '.join(missing)}")

    # Heavy chains only; keep other columns, as real alignments may be among them.
    clonotypes = clonotypes[[c for c in clonotypes.columns if not c.endswith("_light")]]

    # Every exit must write the method file; the workflow fails if it is missing.
    def record_method(method: str, reason: str) -> None:
        if args.out_method is not None:
            args.out_method.write_text(json.dumps({"method": method, "reason": reason}))

    def no_clonotypes() -> None:
        pd.DataFrame(columns=["sequence_id", "clone_id"]).to_csv(
            args.out_clones, sep="\t", index=False
        )
        record_method("none", "no clonotypes for this donor")
        print("no clonotypes for this donor", file=sys.stderr)
        progress("Clustering: no clonotypes")

    # HILARy fails on an empty input, and an empty donor is normal.
    if clonotypes.empty:
        no_clonotypes()
        return

    # No V, J or junction: unplaceable, and HILARy would drop it silently.
    clonotypes = clonotypes[(clonotypes[["v_call", "j_call", "junction"]] != "").all(axis=1)]
    # HILARy trims 3nt per end, so 6nt or less is empty. Renumbered so alignments line up.
    clonotypes = clonotypes[clonotypes["junction"].str.len() > 6].reset_index(drop=True)
    if clonotypes.empty:
        no_clonotypes()
        return

    # HILARy groups by V gene, J gene and CDR3 length, and the adaptive methods raise a
    # KeyError when no group has two members. Each clonotype is then its own lineage.
    if args.mode == "adaptive":
        gene = lambda calls: calls.str.split(",", n=1).str[0].str.split("*", n=1).str[0]
        largest = clonotypes.groupby(
            [gene(clonotypes["v_call"]), gene(clonotypes["j_call"]), clonotypes["junction"].str.len()]
        ).size().max()
        if largest < 2:
            reason = "no lineage has more than one clonotype"
            record_method("singletons", reason)
            print(f"{reason}: no V gene, J gene and CDR3 length group holds two clonotypes",
                  file=sys.stderr)
            progress("Clustering: every clonotype is its own lineage")
            pd.DataFrame({
                "sequence_id": clonotypes["sequence_id"],
                "clone_id": [f"{args.clone_prefix}{i + 1}" for i in range(len(clonotypes))],
            }).to_csv(args.out_clones, sep="\t", index=False)
            return

    work_dir = args.out_clones.parent / "hilary"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Adaptive: full method if every clonotype has an alignment, else cdr3 method.
    method = "crude"
    reason = ""
    alt = None
    if args.mode == "adaptive":
        aligned = (
            pd.read_csv(args.aligned, sep="\t", dtype=str, keep_default_na=False)
            if args.aligned is not None and args.aligned.exists() else pd.DataFrame()
        )
        alt = _alt_alignments(clonotypes, aligned)
        if alt is not None:
            method = "full"
        else:
            method = "cdr3"
            have = clonotypes["sequence_id"].isin(
                set(aligned["sequence_id"]) if not aligned.empty else set())
            reason = (f"{int((~have).sum())} of {len(clonotypes)} clonotypes have no heavy chain "
                      f"alignment, so the shared-mutation test is skipped")
            print(reason, file=sys.stderr)
    record_method(method, reason)

    # No percentage: HILARy reports no usable progress.
    progress(f"Clustering {len(clonotypes)} clonotypes, HILARy {method} method")

    # HILARy mangles ids (`str.strip("-igh")` in, "-igh" appended out), so it gets
    # "s<row>": not numeric, and safe from that strip. Mapped back after.
    surrogates = clonotypes.copy()
    surrogates["_row"] = "s" + surrogates.index.astype(str)
    real_keys = dict(zip(surrogates["_row"], surrogates["sequence_id"]))

    for_hilary = surrogates.drop(columns=["sequence_id"]).rename(columns={"_row": "sequence_id"})
    if alt is not None:
        for_hilary = for_hilary.join(alt.reset_index(drop=True))
    hilary_input = work_dir / "clonotypes.tsv"
    _write_hilary_input(for_hilary, hilary_input)

    clustered = pd.read_csv(
        run_hilary(hilary_input, work_dir, method, args.threshold, args.precision,
                   args.sensitivity, args.threads),
        sep="\t",
    )
    if HILARY_CLONE_COLUMN not in clustered.columns:
        raise RuntimeError(f"HILARy output has no '{HILARY_CLONE_COLUMN}' column")

    returned = clustered["sequence_id"].astype(str).str.removesuffix("-igh")
    unknown = set(returned) - set(real_keys)
    if unknown:
        raise RuntimeError(f"HILARy returned unrecognised sequence ids: {sorted(unknown)[:5]}")

    # String ids, prefixed by donor: HILARy numbers from 1 in every run.
    clones = pd.DataFrame(
        {
            "sequence_id": returned.map(real_keys),
            "clone_id": args.clone_prefix + clustered[HILARY_CLONE_COLUMN].astype(str),
        }
    )
    clones.to_csv(args.out_clones, sep="\t", index=False)


def _read_donor_file(directory: Path, index: int, columns: list[str]) -> pd.DataFrame:
    """Read one donor's table, tolerating a donor that produced nothing."""
    path = directory / f"donor-{index}.tsv"
    if not path.exists():
        return pd.DataFrame(columns=columns)
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return frame if not frame.empty else pd.DataFrame(columns=columns)


def _add_member_nodes(lineages: pd.DataFrame, nodes: pd.DataFrame, links: pd.DataFrame,
                      builders: pd.DataFrame, clonotypes: pd.DataFrame
                      ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Give each lineage without a tree one parentless node per distinct sequence, carrying
    why it has no tree, so its members can still be listed. Tree nodes get a blank reason."""
    nodes = nodes.assign(no_tree_reason="")
    members = lineages[~lineages["lineage_id"].isin(set(nodes["lineage_id"]))]
    if members.empty:
        return nodes, links
    if "group_id" not in members.columns:
        members = members.assign(group_id=members["sequence_id"])
    known = (set(clonotypes.loc[clonotypes["is_known"] == "true", "sequence_id"])
               if "is_known" in clonotypes.columns else set())
    # As the tree step picks a tip's clonotype: a known antibody first, then the lowest id.
    members = (members[["lineage_id", "group_id", "sequence_id"]]
               .drop_duplicates(["lineage_id", "sequence_id"])
               .assign(_other=lambda d: ~d["sequence_id"].isin(known))
               .sort_values(["lineage_id", "group_id", "_other", "sequence_id"]))
    members["is_representative"] = ~members.duplicated(["lineage_id", "group_id"])
    heads = members[members["is_representative"]].copy()
    heads["node_id"] = (heads.groupby("lineage_id").cumcount() + 1).astype(str)
    reason = (builders.set_index("lineage_id")["no_tree_reason"] if not builders.empty
              else pd.Series(dtype=str))
    added = pd.DataFrame({column: "" for column in NODE_COLUMNS}, index=heads.index)
    added["lineage_id"] = heads["lineage_id"]
    added["node_id"] = heads["node_id"]
    added["is_observed"] = "true"
    added["label"] = heads["sequence_id"]
    added["no_tree_reason"] = (heads["lineage_id"].map(reason).replace("", pd.NA)
                               .fillna("No tree was built; the trees log has the details."))
    linked = members.merge(heads[["lineage_id", "group_id", "node_id"]], on=["lineage_id", "group_id"])
    linked = linked.assign(link="1", is_representative=linked["is_representative"].map(
        {True: "true", False: "false"}))[NODE_LINK_FILE_COLUMNS]
    return (pd.concat([nodes, added], ignore_index=True),
            pd.concat([links, linked], ignore_index=True))


def _concat(frames: list[pd.DataFrame], columns: list[str]) -> pd.DataFrame:
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True)


# The membership linker, as exported.
LINEAGE_COLUMNS = ["sequence_id", "lineage_id", "link"]
# As the tree step writes it; group_id sizes a lineage and is not exported.
LINEAGE_FILE_COLUMNS = [*LINEAGE_COLUMNS, "group_id"]
# Branch changes per chain and alphabet.
STEP_COLUMNS = [f"{chain}_{what}" for chain in ("heavy", "light") for what in (
    "mutations_from_parent", "mutation_count_from_parent", "unresolved_from_parent",
    "aa_mutations_from_parent", "aa_mutation_count_from_parent")]
# Changes since the germline and the MRCA, per chain and alphabet, and heavy V and J identity.
HISTORY_COLUMNS = ["distance_from_germline", *[f"{chain}_{what}" for chain in ("heavy", "light") for what in (
    "mutations_from_germline", "mutation_count_from_germline", "mutation_rate_from_germline",
    "aa_mutations_from_germline", "aa_mutation_count_from_germline", "aa_mutation_rate_from_germline",
    "mutations_from_mrca", "mutation_count_from_mrca",
    "aa_mutations_from_mrca", "aa_mutation_count_from_mrca")], "heavy_v_identity", "heavy_j_identity"]
NODE_COLUMNS = ["lineage_id", "node_id", "parent_id", "distance", "is_observed", "label",
                "heavy_sequence", "light_sequence", "node_depth",
                "terminal_branch_fraction", "parent_descendant_count", *STEP_COLUMNS, *HISTORY_COLUMNS]
# The node linker, as exported.
NODE_LINK_COLUMNS = ["lineage_id", "node_id", "sequence_id", "link"]
# As the tree step writes it; is_representative is read here, not exported.
NODE_LINK_FILE_COLUMNS = [*NODE_LINK_COLUMNS, "is_representative"]
BUILDER_COLUMNS = ["lineage_id", "tree_builder", "no_tree_reason"]
GERMLINE_MUTATION_COLUMNS = ["sequence_id", "lineage_id", "germline_mutation_count"]
GERMLINE_MUTATION_OUT_COLUMNS = ["sequence_id", "germline_mutation_count"]
# Tooltip content per observed node. Always all written, blank if absent.
NODE_PROPERTY_SOURCE_COLUMNS = ["v_call", "j_call", "junction", "cdr1_aa", "cdr2_aa", "cdr3_aa",
                                "sequence_aa", "main_sequence",
                                "v_call_light", "j_call_light", "junction_light", "cdr1_aa_light",
                                "cdr2_aa_light", "cdr3_aa_light", "sequence_aa_light",
                                "main_sequence_light", "is_known"]
KNOWN_DISTANCE_COLUMNS = ["sequence_id", "known_id", "known_aa_heavy", "known_nt_heavy",
                           "known_aa_light", "known_nt_light"]
NODE_PROPERTY_COLUMNS = ["lineage_id", "node_id", *NODE_PROPERTY_SOURCE_COLUMNS]
# Every count a dataset gives beside its primary abundance, by the abundance table's header.
# Blank where no dataset of a node or lineage carries that count.
COUNT_COLUMNS = ["reads", "umis", "cells"]


def _count_sums(frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Each count summed per key, blank where no row has it; whole numbers stay whole."""
    counts = [c for c in COUNT_COLUMNS if c in frame.columns]
    summed = frame.groupby(keys, as_index=False)[counts].sum(min_count=1) if counts else (
        frame[keys].drop_duplicates())
    for column in COUNT_COLUMNS:
        if column not in summed.columns:
            summed[column] = float("nan")
        values = summed[column]
        whole = values.dropna()
        if (whole == whole.round(0)).all():
            summed[column] = values.round(0).astype("Int64")
    return summed[[*keys, *COUNT_COLUMNS]]



def _one_per_clonotype(frame: pd.DataFrame, sizes: pd.Series) -> pd.DataFrame:
    """One row per clonotype, so lead selection can rank by it.

    A clonotype in two donors' lineages keeps the one from the larger lineage.
    """
    if frame.empty:
        return frame
    ranked = frame.assign(_support=frame["lineage_id"].map(sizes).fillna(0).astype(int))
    ranked = ranked.sort_values(["sequence_id", "_support", "lineage_id"],
                                ascending=[True, False, True])
    deduped = ranked.drop_duplicates("sequence_id", keep="first")
    if len(deduped) < len(ranked):
        print(f"{len(ranked) - len(deduped)} clonotypes were scored in more than one "
              f"donor's lineage; kept the better-supported one")
    return deduped.drop(columns="_support")


def _write_one_per_clonotype(frame: pd.DataFrame, sizes: pd.Series, path: Path,
                             columns: list[str]) -> None:
    """Write one row per clonotype, no sample axis."""
    _one_per_clonotype(frame, sizes).reindex(columns=columns).to_csv(path, sep="\t", index=False)


def _gene(call: str) -> str:
    """Gene-level name of a call: first hit, allele stripped."""
    return call.split(",")[0].split("*")[0]


def _mode(values: pd.Series) -> str:
    """Most common non-empty value, ties broken alphabetically; blank if none."""
    present = values[values.astype(str).str.len() > 0]
    if present.empty:
        return ""
    counts = present.value_counts()
    return sorted(counts[counts == counts.max()].index)[0]


def _lineage_descriptors(
    lineages: pd.DataFrame,
    clonotypes: pd.DataFrame,
    present: pd.DataFrame | None,
    dataset_total: float | None,
) -> pd.DataFrame:
    """Per lineage: modal genes and CDR3 length, the top member's CDR3, and abundance.

    The abundance fraction is of the whole dataset.
    """
    members = lineages.merge(clonotypes, on="sequence_id", how="left")
    for column in ("v_call", "j_call", "junction", "cdr3_aa", "v_call_light", "j_call_light"):
        if column not in members.columns:
            members[column] = ""
        members[column] = members[column].fillna("").astype(str)
    members["v_gene"] = members["v_call"].map(_gene)
    members["j_gene"] = members["j_call"].map(_gene)
    members["light_v_gene"] = members["v_call_light"].map(_gene)
    members["light_j_gene"] = members["j_call_light"].map(_gene)
    # MiXCR's CDR3 length, Cys and Trp included.
    members["cdr3_length_aa"] = members["junction"].str.len() // 3

    if present is not None and not present.empty:
        weight = present.groupby(["lineage_id", "sequence_id"])["abundance"].sum()
        members["weight"] = [
            weight.get((lineage, sequence), 0.0)
            for lineage, sequence in zip(members["lineage_id"], members["sequence_id"])
        ]
    else:
        members["weight"] = 0.0
    members = members.sort_values(["lineage_id", "weight", "sequence_id"],
                                  ascending=[True, False, True], kind="stable")

    rows = []
    for lineage_id, group in members.groupby("lineage_id", sort=True):
        rows.append({
            "lineage_id": lineage_id,
            "v_gene": _mode(group["v_gene"]),
            "j_gene": _mode(group["j_gene"]),
            "light_v_gene": _mode(group["light_v_gene"]),
            "light_j_gene": _mode(group["light_j_gene"]),
            "cdr3_length_aa": _mode(group["cdr3_length_aa"].astype(str)),
            "representative_cdr3_aa": group["cdr3_aa"].iloc[0],
        })
    described = pd.DataFrame(rows, columns=[
        "lineage_id", "v_gene", "j_gene", "light_v_gene", "light_j_gene",
        "cdr3_length_aa", "representative_cdr3_aa",
    ])

    if present is not None and not present.empty:
        totals = present.groupby("lineage_id").agg(
            total_abundance=("abundance", "sum"), sample_count=("sample_id", "nunique"),
        )
        described["total_abundance"] = described["lineage_id"].map(totals["total_abundance"]).fillna(0)
        described["sample_count"] = described["lineage_id"].map(totals["sample_count"]).fillna(0).astype(int)
        described["abundance_fraction"] = (
            described["total_abundance"] / dataset_total if dataset_total else float("nan")
        )
        whole = described["total_abundance"].round(0)
        described["total_abundance"] = whole.astype(int) if (whole == described["total_abundance"]).all() else described["total_abundance"]
        totals = _count_sums(present, ["lineage_id"]).set_index("lineage_id")
        for column in COUNT_COLUMNS:
            described[f"total_{column}"] = described["lineage_id"].map(totals[column])
    else:
        described["total_abundance"] = ""
        described["sample_count"] = ""
        described["abundance_fraction"] = ""
        for column in COUNT_COLUMNS:
            described[f"total_{column}"] = ""
    return described


def _node_properties(links: pd.DataFrame, clonotypes: pd.DataFrame,
                     lineages: pd.DataFrame) -> pd.DataFrame:
    """Copy each observed node's representative clonotype onto the node's axes.

    The dendrogram cannot follow a linker, so tooltip content must sit on the tree axes.
    Per donor where the tables carry one: alleles are inferred per donor, so a clonotype
    in two donors may have a different V call in each.
    """
    available = [c for c in NODE_PROPERTY_SOURCE_COLUMNS if c in clonotypes.columns]
    if "is_representative" in links.columns:
        links = links[links["is_representative"] == "true"]
    keys = ["sequence_id"]
    links = links[["lineage_id", "node_id", "sequence_id"]]
    if "donor" in clonotypes.columns:
        keys.append("donor")
        donor_of = (lineages.drop_duplicates("lineage_id").set_index("lineage_id")["donor"]
                    if "donor" in lineages.columns else pd.Series(dtype=str))
        links = links.assign(donor=links["lineage_id"].map(donor_of).fillna(""))
    joined = links.merge(clonotypes[[*keys, *available]], on=keys, how="inner")
    for column in NODE_PROPERTY_SOURCE_COLUMNS:
        if column not in joined.columns:
            joined[column] = ""
    return joined[NODE_PROPERTY_COLUMNS]


def _short_labels(nodes: pd.DataFrame, links: pd.DataFrame, clonotypes: pd.DataFrame) -> pd.Series:
    """Node labels, using the representative's short MiXCR clone label where there is one."""
    if "clone_label" not in clonotypes.columns or links.empty:
        return nodes["label"]
    if "is_representative" in links.columns:
        links = links[links["is_representative"] == "true"]
    short = (
        links[["lineage_id", "node_id", "sequence_id"]]
        .merge(clonotypes[["sequence_id", "clone_label"]], on="sequence_id", how="inner")
        .drop_duplicates(["lineage_id", "node_id"])
    )
    short = short[short["clone_label"].astype(str) != ""]
    by_node = short.set_index(["lineage_id", "node_id"])["clone_label"]
    keys = pd.MultiIndex.from_arrays([nodes["lineage_id"].astype(str), nodes["node_id"].astype(str)])
    found = pd.Series(by_node.reindex(keys).to_numpy(), index=nodes.index)
    return found.where(found.notna(), nodes["label"])


def _data_sources(lineages: pd.DataFrame, clonotypes: pd.DataFrame) -> pd.Series:
    """Sorted dataset names of each lineage's members, keyed by lineage id."""
    if "data_source" not in clonotypes.columns or lineages.empty:
        return pd.Series(dtype=str)
    members = lineages[["lineage_id", "sequence_id"]].merge(
        clonotypes[["sequence_id", "data_source"]], on="sequence_id", how="left")
    return members.groupby("lineage_id")["data_source"].agg(
        lambda values: ", ".join(sorted({v for v in values if isinstance(v, str) and v})))


def _write_per_dataset(
    args: argparse.Namespace,
    lineages: pd.DataFrame,
    links: pd.DataFrame,
    known_distances: pd.DataFrame,
    germline_mutations: pd.DataFrame,
    known_labels: pd.Series,
) -> None:
    """Write each dataset's rows with the id prefix removed; every file exists, even if empty."""
    args.per_dataset_dir.mkdir(parents=True, exist_ok=True)
    sizes = (lineages.groupby("lineage_id").size() if not lineages.empty
             else pd.Series(dtype=int))

    def slice_of(frame: pd.DataFrame, index: str, column: str) -> pd.DataFrame:
        if frame.empty:
            return frame.copy()
        part = frame[dataset_of(frame[column]) == index].copy()
        part[column] = unprefixed(part[column])
        return part

    all_lineages = sorted(lineages["lineage_id"].drop_duplicates()) if not lineages.empty else []
    for index in args.dataset:
        out = args.per_dataset_dir
        members = slice_of(lineages, index, "sequence_id")[LINEAGE_COLUMNS]
        members.to_csv(out / f"lineages-{index}.tsv", sep="\t", index=False)
        # This dataset's members per lineage; every lineage gets a row, zero if none.
        counts = members.groupby("lineage_id").size()
        pd.DataFrame({
            "lineage_id": all_lineages,
            "member_count": pd.Series(all_lineages).map(counts).fillna(0).astype(int).to_numpy(),
        }).to_csv(out / f"member-counts-{index}.tsv", sep="\t", index=False)
        slice_of(links, index, "sequence_id")[NODE_LINK_COLUMNS].to_csv(
            out / f"node-links-{index}.tsv", sep="\t", index=False)
        _write_one_per_clonotype(slice_of(germline_mutations, index, "sequence_id"), sizes,
                                 out / f"germline-mutations-{index}.tsv",
                                 GERMLINE_MUTATION_OUT_COLUMNS)
        distances = _one_per_clonotype(slice_of(known_distances, index, "sequence_id"), sizes)
        # Candidate-to-known antibody links, one file per known antibody dataset, written even if empty.
        for known_index in args.known:
            pair = (
                distances[dataset_of(distances["known_id"]) == known_index]
                if not distances.empty else distances
            )
            linked = pd.DataFrame({
                "sequence_id": pair["sequence_id"] if not pair.empty else pd.Series(dtype=str),
                "known_id": unprefixed(pair["known_id"]) if not pair.empty else pd.Series(dtype=str),
                "link": 1,
            })
            linked.to_csv(
                out / f"known-links-{index}-{known_index}.tsv", sep="\t", index=False)
        # The known antibody's readable label where its producer gives one (MiXCR's Clone Id), else its id
        # as its own dataset knows it. The links above keep the id: it is their join key.
        if not distances.empty:
            label = distances["known_id"].map(known_labels)
            distances["known_id"] = label.where(label.notna() & (label != ""),
                                                 unprefixed(distances["known_id"]))
        distances.reindex(columns=KNOWN_DISTANCE_COLUMNS).to_csv(
            out / f"known-distances-{index}.tsv", sep="\t", index=False)


def _read_clonotypes(args: argparse.Namespace, donors: list) -> pd.DataFrame:
    """The clonotypes the tools saw, with merge's dataset name and known antibody flag joined back.

    From per-donor files, one row per donor and clonotype, with a `donor` column.
    """
    clonotypes = _read_clonotype_tables(args, donors)
    annotations = getattr(args, "annotations", None)
    if annotations is None or not annotations.exists():
        return clonotypes
    marks = pd.read_csv(annotations, sep="\t", dtype=str, keep_default_na=False)
    kept = [c for c in clonotypes.columns if c not in ("data_source", "is_known")]
    return clonotypes[kept].merge(marks.drop_duplicates("sequence_id"), on="sequence_id", how="left").fillna("")


def _read_clonotype_tables(args: argparse.Namespace, donors: list) -> pd.DataFrame:
    """Per-donor files (with realigned calls) if given, else the merged file."""
    if args.clonotypes_dir is not None:
        parts = [_read_donor_file(args.clonotypes_dir, index, ["sequence_id"])
                 .assign(donor="" if donor is None else donor)
                 for index, donor in enumerate(donors)]
        frames = [f.fillna("") for f in parts if not f.empty]
        if not frames:
            return pd.DataFrame(columns=["sequence_id", "donor"])
        return pd.concat(frames, ignore_index=True)
    if args.clonotypes is not None:
        return pd.read_csv(args.clonotypes, sep="\t", dtype=str, keep_default_na=False)
    return pd.DataFrame(columns=["sequence_id"])


def _node_abundance(present, links: pd.DataFrame, clonotypes: pd.DataFrame) -> pd.DataFrame:
    """Each node's abundance over all samples and clonotypes, its clonotype count and datasets.

    No sample axis: the dendrogram could not join it.
    """
    columns = ["lineage_id", "node_id", "abundance", "clonotype_count", "dataset", *COUNT_COLUMNS]
    if present is None or links.empty:
        return pd.DataFrame(columns=columns)
    counts = [c for c in COUNT_COLUMNS if c in present.columns]
    placed = links[["lineage_id", "node_id", "sequence_id"]].merge(
        present[["sequence_id", "lineage_id", "abundance", *counts]], on=["lineage_id", "sequence_id"],
    )
    if placed.empty:
        return pd.DataFrame(columns=columns)
    summed = placed.groupby(["lineage_id", "node_id"], as_index=False)["abundance"].sum().merge(
        _count_sums(placed, ["lineage_id", "node_id"]), on=["lineage_id", "node_id"])
    # Datasets as merge named them, "A, B" when seen in both. One flag per dataset rides the
    # count's groupby; single-dataset runs skip it.
    names = []
    if "data_source" in clonotypes.columns:
        names = sorted(set(clonotypes["data_source"]) - {""})
    flags = {}
    if len(names) > 1:
        source = clonotypes.drop_duplicates("sequence_id").set_index("sequence_id")["data_source"]
        of_link = source.reindex(links["sequence_id"]).to_numpy()
        flags = {f"_in{i}": of_link == name for i, name in enumerate(names)}
    counted = (
        links[["lineage_id", "node_id", "sequence_id"]].assign(**flags)
        .groupby(["lineage_id", "node_id"], as_index=False)
        .agg(clonotype_count=("sequence_id", "nunique"), **{f: (f, "max") for f in flags})
    )
    if flags:
        mask = sum(counted[f"_in{i}"].astype("int64") * (1 << i) for i in range(len(names)))
        labels = {m: ", ".join(n for i, n in enumerate(names) if m >> i & 1) for m in mask.unique()}
        counted = counted.drop(columns=list(flags)).assign(dataset=mask.map(labels))
    return summed.merge(counted, on=["lineage_id", "node_id"], how="right").fillna(
        {"abundance": 0, "dataset": ""},
    ).reindex(columns=columns)


def _isotypes(clonotypes: pd.DataFrame) -> pd.Series:
    """Each clonotype's isotype by id: the producer's, else its heavy C gene's class (IGHG3 is IgG)."""
    given = (clonotypes["isotype"].fillna("").astype(str) if "isotype" in clonotypes.columns
             else pd.Series("", index=clonotypes.index))
    if "c_call" in clonotypes.columns:
        cls = clonotypes["c_call"].fillna("").astype(str).str.extract(r"^IGH([ADEGM])", expand=False)
        given = given.where(given != "", ("Ig" + cls).fillna(""))
    return pd.Series(given.to_numpy(), index=clonotypes["sequence_id"].astype(str))


# Counts an isotype vote is weighed in, most deduplicated first.
ISOTYPE_WEIGHTS = ["cells", "umis", "reads"]


def _node_isotype(present, links: pd.DataFrame, clonotypes: pd.DataFrame) -> pd.DataFrame:
    """Each observed node's isotype, ties alphabetical.

    A node is one sequence, and clonotypes split by C gene share it. One isotype among them
    wins outright. Otherwise the isotypes are weighed in the first count every clonotype of the
    node has (datasets count in different units, which cannot be compared); with none, the node
    is left blank. Without abundance each clonotype weighs one. Nodes without any isotype are
    left out.
    """
    columns = ["lineage_id", "node_id", "isotype"]
    isotype = _isotypes(clonotypes.drop_duplicates("sequence_id"))
    node = ["lineage_id", "node_id"]
    placed = links[[*node, "sequence_id"]].drop_duplicates()
    placed = placed.assign(isotype=placed["sequence_id"].map(isotype).fillna(""))
    placed = placed[placed["isotype"] != ""].reset_index(drop=True)
    if placed.empty:
        return pd.DataFrame(columns=columns)
    single = placed.groupby(node)["isotype"].transform("nunique") == 1
    if present is None or present.empty:
        placed["weight"] = 1.0
    else:
        units = [u for u in ISOTYPE_WEIGHTS if u in present.columns]
        placed["weight"] = float("nan")
        if units:
            per = present.groupby(["lineage_id", "sequence_id"])[units].sum(min_count=1)
            keys = pd.MultiIndex.from_arrays([placed["lineage_id"], placed["sequence_id"]])
            values = per.reindex(keys).reset_index(drop=True)
            open_ = placed["weight"].isna()
            for unit in units:
                shared = values[unit].notna().groupby([placed["lineage_id"], placed["node_id"]]).transform("all")
                take = open_ & shared
                placed.loc[take, "weight"] = values.loc[take, unit]
                open_ &= ~shared
        placed.loc[single, "weight"] = placed.loc[single, "weight"].fillna(1.0)
    placed = placed[placed["weight"].notna()]
    summed = placed.groupby([*node, "isotype"], as_index=False)["weight"].sum()
    summed = summed.sort_values([*node, "weight", "isotype"],
                                ascending=[True, True, False, True], kind="stable")
    return summed.drop_duplicates(node)[columns]


def _node_metadata(args: argparse.Namespace, present, links: pd.DataFrame) -> pd.DataFrame:
    """Each node's distinct metadata values joined with ", ", plus a <column>__count of them."""
    parts = []
    for index, path in args.sample_metadata:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        table["sample_id"] = prefixed(index, table["sample_id"])
        parts.append(table)
    meta_columns = list(dict.fromkeys(c for t in parts for c in t.columns if c != "sample_id"))
    columns = ["lineage_id", "node_id"] + meta_columns
    # Every exit writes the same header: the workflow imports each count column.
    header = columns + [f"{c}__count" for c in meta_columns]
    if present is None or links.empty or not parts:
        return pd.DataFrame(columns=header)
    metadata = pd.concat(parts, ignore_index=True).fillna("")
    placed = (
        present.loc[present["abundance"] > 0, ["sample_id", "sequence_id", "lineage_id"]]
        .merge(links[["lineage_id", "sequence_id", "node_id"]], on=["lineage_id", "sequence_id"])
        .merge(metadata, on="sample_id")
    )
    if placed.empty:
        return pd.DataFrame(columns=header)
    distinct = lambda values: ", ".join(dict.fromkeys(v for v in values if v))
    counted = lambda values: len(dict.fromkeys(v for v in values if v))
    merged = (
        placed.groupby(["lineage_id", "node_id"], as_index=False)
        .agg({column: distinct for column in meta_columns})
    )
    counts = (
        placed.groupby(["lineage_id", "node_id"], as_index=False)
        .agg({column: counted for column in meta_columns})
        .rename(columns={column: f"{column}__count" for column in meta_columns})
    )
    return merged.merge(counts, on=["lineage_id", "node_id"]).reindex(columns=header)


def _run_id(lineages: pd.DataFrame, nodes: pd.DataFrame) -> str:
    """Content id of lineage membership and tree shapes, so saved views detect renumbering."""
    digest = hashlib.sha256()
    for frame, columns in ((lineages, ["sequence_id", "lineage_id"]),
                           (nodes, ["lineage_id", "node_id", "parent_id", "label"])):
        part = frame.reindex(columns=columns).fillna("").astype(str).sort_values(columns)
        digest.update(part.to_csv(sep="\t", index=False).encode())
    return digest.hexdigest()[:16]


def collect(args: argparse.Namespace) -> None:
    """Merge every donor's tree output, read by position as `split` wrote it, and summarise."""
    # No donor column: one unnamed group at position 0.
    donors = args.donor or [None]
    steps = ["Reading tree outputs", "Labelling nodes", "Summarising lineages",
             "Abundance and metadata", "Describing lineages", "Writing per-dataset tables"]

    def step(text: str) -> None:
        progress(f"{text}: {100 * steps.index(text) / len(steps):.0f}%")

    step("Reading tree outputs")

    lineage_parts, node_parts, link_parts, builder_parts = [], [], [], []
    distance_parts, germline_parts = [], []
    for index, donor in enumerate(donors):
        lineages = _read_donor_file(args.lineages_dir, index, LINEAGE_FILE_COLUMNS)
        if donor is not None:
            lineages = lineages.assign(donor=donor)
        lineage_parts.append(lineages)
        node_parts.append(_read_donor_file(args.nodes_dir, index, NODE_COLUMNS))
        link_parts.append(_read_donor_file(args.node_links_dir, index, NODE_LINK_FILE_COLUMNS))
        if args.builders_dir is not None:
            builder_parts.append(_read_donor_file(args.builders_dir, index, BUILDER_COLUMNS))
        if args.germline_mutations_dir is not None:
            germline_parts.append(_read_donor_file(args.germline_mutations_dir, index,
                                                   GERMLINE_MUTATION_COLUMNS))
        if args.known_distances_dir is not None:
            # With the donor's lineage, so a clonotype scored in two donors keeps one score.
            distances = _read_donor_file(args.known_distances_dir, index, KNOWN_DISTANCE_COLUMNS)
            distance_parts.append(distances.merge(
                lineages[["sequence_id", "lineage_id"]].drop_duplicates("sequence_id"),
                on="sequence_id", how="left"))

    lineages = _concat(lineage_parts, LINEAGE_FILE_COLUMNS)
    # Empty donors leave no lineage rows, and the concat then drops the donor column too.
    if args.donor and "donor" not in lineages.columns:
        lineages["donor"] = pd.Series(dtype=str)
    nodes = _concat(node_parts, NODE_COLUMNS)
    # Tip labels are prefixed clonotype ids; strip the prefix.
    links = _concat(link_parts, NODE_LINK_FILE_COLUMNS)
    donor_clonotypes = _read_clonotypes(args, donors)
    # One row per clonotype for everything but the node properties, which keep each donor's calls.
    clonotypes = donor_clonotypes.drop_duplicates("sequence_id").drop(columns="donor", errors="ignore")
    builders = _concat(builder_parts, BUILDER_COLUMNS)
    nodes, links = _add_member_nodes(lineages, nodes, links, builders, clonotypes)
    step("Labelling nodes")
    if not nodes.empty:
        observed_label = nodes["label"].astype(str).str.contains(DATASET_SEP, regex=False)
        nodes.loc[observed_label, "label"] = unprefixed(nodes.loc[observed_label, "label"])
        nodes["label"] = _short_labels(nodes, links, clonotypes)
    nodes.to_csv(args.out_nodes, sep="\t", index=False)
    if args.out_node_properties is not None:
        _node_properties(links, donor_clonotypes, lineages).to_csv(
            args.out_node_properties, sep="\t", index=False)

    step("Summarising lineages")
    # cluster_size: distinct sequences (tree step groups). tip_count: tips actually
    # drawn, zero when no tree was built. Without groups, each clonotype is its own.
    sized = lineages.copy()
    if "group_id" not in sized.columns:
        sized["group_id"] = sized["sequence_id"]
    lineage_stats = (
        sized.groupby("lineage_id", as_index=False)["group_id"]
        .nunique()
        .rename(columns={"group_id": "cluster_size"})
    )
    # Tree tips only: the member nodes of lineages without a tree carry a reason.
    observed = (nodes[(nodes["is_observed"] == "true") & (nodes["no_tree_reason"] == "")]
                if not nodes.empty else nodes)
    tips = observed.groupby("lineage_id").size() if not observed.empty else pd.Series(dtype=int)
    lineage_stats["tip_count"] = (
        lineage_stats["lineage_id"].map(tips).fillna(0).astype(int)
    )

    # The builder sets the branch length unit: steps or substitutions per site.
    by_lineage = builders.set_index("lineage_id") if not builders.empty else None
    for column in ("tree_builder", "no_tree_reason"):
        lineage_stats[column] = lineage_stats["lineage_id"].map(
            by_lineage[column] if by_lineage is not None else pd.Series(dtype=str))

    # TODO(badges): the pass / alert / ignore quality badge (spec Deliverable 1) goes here.

    step("Abundance and metadata")
    present = None
    dataset_total = None
    have_abundance = args.abundance is not None and args.abundance.exists()
    if have_abundance and not lineages.empty:
        abundance = pd.read_csv(
            args.abundance,
            sep="\t",
            dtype={"sample_id": str, "sequence_id": str},
            keep_default_na=False,
        )
        if "donor" in lineages.columns:
            # Join on donor too, so a clonotype in two donors credits each lineage its own counts.
            sample_donors = pd.read_csv(args.donors, sep="\t", dtype=str, keep_default_na=False)
            abundance = abundance.merge(sample_donors, on="sample_id", how="inner")
            present = abundance.merge(lineages, on=["sequence_id", "donor"], how="inner")
        else:
            present = abundance.merge(lineages, on="sequence_id", how="inner")
        dataset_total = float(pd.to_numeric(abundance["abundance"], errors="coerce").fillna(0).sum())
        present = present.assign(abundance=pd.to_numeric(present["abundance"], errors="coerce").fillna(0))
        # Blank stays blank: a dataset without the count must not read as zero.
        for column in COUNT_COLUMNS:
            if column in present.columns:
                present[column] = pd.to_numeric(present[column], errors="coerce")
    if args.out_node_metadata is not None:
        _node_metadata(args, present, links).to_csv(args.out_node_metadata, sep="\t", index=False)
    if args.out_node_abundance is not None:
        abundance = _node_abundance(present, links, clonotypes).merge(
            _node_isotype(present, links, clonotypes), on=["lineage_id", "node_id"], how="left")
        abundance.fillna({"isotype": ""}).to_csv(args.out_node_abundance, sep="\t", index=False)
    step("Describing lineages")
    described = _lineage_descriptors(lineages, clonotypes, present, dataset_total)
    lineage_stats = lineage_stats.merge(described, on="lineage_id", how="left")
    # Show the donor in its own column and the label without it; the id is unchanged.
    if "donor" in lineages.columns:
        donor_of = lineages.drop_duplicates("lineage_id").set_index("lineage_id")["donor"]
        lineage_stats["donor"] = lineage_stats["lineage_id"].map(donor_of).fillna("")
        lineage_stats["lineage_label"] = [
            lineage[len(donor) + 1:] if donor and lineage.startswith(donor + "/") else lineage
            for lineage, donor in zip(lineage_stats["lineage_id"], lineage_stats["donor"])
        ]
    # The full id, exported as the lineage's label so other blocks can show it.
    lineage_stats["lineage_name"] = lineage_stats["lineage_id"]
    lineage_stats["data_source"] = (
        lineage_stats["lineage_id"].map(_data_sources(lineages, clonotypes)).fillna("")
    )
    # Known antibodies per lineage, for the known antibody search.
    if "is_known" in clonotypes.columns and not lineages.empty:
        known = set(clonotypes.loc[clonotypes["is_known"] == "true", "sequence_id"])
        held = lineages[lineages["sequence_id"].isin(known)].groupby("lineage_id").size()
        lineage_stats["known_count"] = lineage_stats["lineage_id"].map(held).fillna(0).astype(int)
    else:
        lineage_stats["known_count"] = 0
    lineage_stats.to_csv(args.out_lineage_stats, sep="\t", index=False)

    step("Writing per-dataset tables")
    known_labels = (
        clonotypes.drop_duplicates("sequence_id").set_index("sequence_id")["clone_label"]
        if "clone_label" in clonotypes.columns else pd.Series(dtype=str))
    _write_per_dataset(args, lineages, links,
                       _concat(distance_parts, [*KNOWN_DISTANCE_COLUMNS, "lineage_id"]),
                       _concat(germline_parts, GERMLINE_MUTATION_COLUMNS), known_labels)

    # Per-donor counts for the overview; names match the workflow's groups.
    if args.out_donor_stats is not None:
        stats = []
        for donor in donors:
            own = lineages[lineages["donor"] == donor] if donor is not None else lineages
            stats.append({
                "donor": donor if donor is not None else "all samples",
                "clonotype_count": int(own["sequence_id"].nunique()) if not own.empty else 0,
                "lineage_count": int(own["lineage_id"].nunique()) if not own.empty else 0,
            })
        args.out_donor_stats.write_text(json.dumps(stats))

    if args.out_run_id is not None:
        args.out_run_id.write_text(json.dumps({"runId": _run_id(lineages, nodes)}))

    if args.logs_dir is not None and args.out_log is not None:
        with args.out_log.open("w") as out:
            for index, donor in enumerate(donors):
                path = args.logs_dir / f"donor-{index}.log"
                out.write(f"=== {donor if donor is not None else 'all samples'} ===\n")
                out.write(path.read_text() if path.exists() else "no log\n")
    progress("Collected: 100%")


def count_hits(args: argparse.Namespace) -> None:
    """Count the nodes and lineages a sequence search hit, as JSON for the UI."""
    hits = pd.read_csv(args.hits, sep="\t", dtype=str, keep_default_na=False)
    lineage_column = hits.columns[0] if len(hits.columns) else None
    args.out.write_text(json.dumps({
        "nodes": int(len(hits)),
        "lineages": int(hits[lineage_column].nunique()) if lineage_column and len(hits) else 0,
    }))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest="stage", required=True)

    m = stages.add_parser("merge", help="every dataset's tables as one, ids prefixed by dataset")
    m.add_argument(
        "--dataset",
        action="append",
        nargs=4,
        metavar=("INDEX", "CLONOTYPES", "ABUNDANCE", "NAME"),
        required=True,
        help="a dataset's position, its clonotype table, its abundance table, "
             "and its name as the lineage table shows it",
    )
    m.add_argument(
        "--known",
        action="append",
        default=[],
        help="position of a dataset whose clonotypes are known antibodies; repeat per such dataset",
    )
    m.add_argument("--out-clonotypes", required=True, type=Path)
    m.add_argument("--out-abundance", required=True, type=Path)
    m.add_argument("--out-annotations", type=Path,
                   help="dataset name and known antibody flag per clonotype, kept out of --out-clonotypes")
    m.add_argument("--out-known", type=Path,
                   help="sequence_id and is_known of the known antibody clonotypes only, for the tree step")
    m.add_argument("--out-gene-usage", type=Path,
                   help="JSON: clonotypes per heavy V and J gene, per dataset")
    m.set_defaults(func=merge)

    sp = stages.add_parser("split", help="one clonotype table per donor")
    sp.add_argument("--clonotypes", required=True, type=Path)
    sp.add_argument("--abundance", type=Path)
    sp.add_argument("--donors", type=Path, help="sample_id, donor")
    sp.add_argument("--out-dir", required=True, type=Path)
    sp.add_argument(
        "--donor",
        action="append",
        default=[],
        help="donor name; repeat once per donor, in the order the output files are numbered",
    )
    sp.set_defaults(func=split)

    c = stages.add_parser("cluster", help="assign heavy-chain clone ids")
    c.add_argument("--clonotypes", required=True, type=Path)
    c.add_argument("--out-clones", required=True, type=Path)
    c.add_argument("--mode", choices=["fixed", "adaptive"], default="fixed",
                   help="fixed: crude method at --threshold; adaptive: full method at "
                        "--precision/--sensitivity, or cdr3 method where alignments are missing")
    c.add_argument(
        "--threshold",
        default=0.2,
        type=float,
        help="fixed mode: fraction of CDR3 length below which two clonotypes share a clone",
    )
    c.add_argument("--precision", default=0.99, type=float, help="adaptive mode: desired precision")
    c.add_argument("--sensitivity", default=0.9, type=float, help="adaptive mode: desired sensitivity")
    c.add_argument("--aligned", type=Path,
                   help="the align stage's table, for the full method's mutation counts")
    c.add_argument("--out-method", type=Path, help="JSON: which HILARy method ran, and why")
    c.add_argument("--threads", default=1, type=int)
    c.add_argument(
        "--clone-prefix",
        default="",
        help="prepended to every clone id, so donors cannot collide on one",
    )
    c.set_defaults(func=cluster)

    h = stages.add_parser("count-hits", help="nodes and lineages a sequence search hit, as JSON")
    h.add_argument("--hits", required=True, type=Path, help="the search's hits table, lineage key first")
    h.add_argument("--out", required=True, type=Path)
    h.set_defaults(func=count_hits)

    k = stages.add_parser("collect", help="merge per-donor tree output and summarise")
    k.add_argument("--lineages-dir", required=True, type=Path)
    k.add_argument("--annotations", type=Path, help="merge's dataset name and known antibody flag per clonotype")
    k.add_argument(
        "--donor",
        action="append",
        default=[],
        help="donor name; repeat in the same order `split` was given",
    )
    k.add_argument("--nodes-dir", required=True, type=Path)
    k.add_argument("--node-links-dir", required=True, type=Path)
    k.add_argument("--builders-dir", type=Path)
    k.add_argument("--germline-mutations-dir", type=Path)
    k.add_argument("--known-distances-dir", type=Path)
    k.add_argument("--abundance", type=Path)
    k.add_argument("--donors", type=Path)
    k.add_argument("--out-nodes", required=True, type=Path)
    k.add_argument("--out-lineage-stats", required=True, type=Path)
    k.add_argument("--clonotypes", type=Path,
                   help="the merged clonotype table, for the node properties, the lineage "
                        "descriptors and the datasets each lineage drew members from")
    k.add_argument("--clonotypes-dir", type=Path,
                   help="instead of --clonotypes: donor-<i>.tsv tables as the tools saw them")
    k.add_argument("--out-node-properties", type=Path)
    k.add_argument("--sample-metadata", action="append", nargs=2, default=[], metavar=("INDEX", "PATH"),
                   help="a dataset's sample_id plus meta_<k> columns; repeat once per dataset")
    k.add_argument("--out-node-metadata", type=Path)
    k.add_argument(
        "--known",
        action="append",
        default=[],
        help="position of a dataset whose clonotypes are known antibodies; repeat per such dataset",
    )
    k.add_argument("--out-node-abundance", type=Path,
                   help="each node's summed abundance and how many clonotypes it stands for")
    k.add_argument("--out-donor-stats", type=Path, help="per-donor clonotype and lineage counts, JSON")
    k.add_argument("--out-run-id", type=Path, help="content id of the lineages and trees, JSON")
    k.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="a dataset position `merge` prefixed ids with; repeat once per dataset",
    )
    k.add_argument("--per-dataset-dir", required=True, type=Path,
                   help="where each dataset's slice of the per-clonotype and per-sample tables goes")
    k.add_argument("--logs-dir", type=Path)
    k.add_argument("--out-log", type=Path)
    k.set_defaults(func=collect)

    args = parser.parse_args()
    args.func(args)
    print(f"{args.stage} finished in {clock()}", flush=True)


if __name__ == "__main__":
    main()
