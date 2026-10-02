#!/usr/bin/env Rscript
# Checks for trees.R. Run through test/run.sh, which builds the fixtures first.
# Light chains refine lineages but never cost a clonotype its lineage; light-less clonotypes join their paired kin.

root <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(root)) stop("usage: test_trees.R <fixture dir>")
trees_R <- "/app/trees.R"

failures <- 0L
checks <- 0L
ok <- function(what, cond) {
  checks <<- checks + 1L
  if (isTRUE(cond)) cat(sprintf("  ok   %s\n", what))
  else { cat(sprintf("  FAIL %s\n", what)); failures <<- failures + 1L }
}

read_tsv <- function(p) read.delim(p, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

# Runs the align and trees stages as the workflow does; both logs are joined.
# `airr = FALSE` passes no AIRR dir, as for imported data.
run_one <- function(scenario, ..., light = FALSE, builder = "raxml", airr = TRUE) {
  dir <- file.path(root, scenario)
  outdir <- file.path(dir, paste0("out-", builder, if (light) "-light" else ""))
  dir.create(outdir, showWarnings = FALSE)
  aligned <- file.path(outdir, "aligned.tsv")
  common <- c("--clonotypes", file.path(dir, "clonotypes.tsv"),
              if (light) "--light-chains")
  align_argv <- c(trees_R, "--stage", "align", common,
                  if (airr) c("--airr-dir", file.path(dir, "airr")),
                  "--out-aligned", aligned)
  align_log <- system2("Rscript", align_argv, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(align_log, "status"))) {
    return(list(ok = FALSE, log = paste(align_log, collapse = "\n"), dir = outdir))
  }
  marks <- file.path(dir, "annotations.tsv")
  argv <- c(trees_R, "--stage", "trees", common,
            if (file.exists(marks)) c("--annotations", marks),
            "--aligned", aligned,
            "--clones", file.path(dir, "clones.tsv"),
            # IgPhyML is chosen by scope; trees.R has no --builder.
            if (builder == "igphyml") c("--igphyml-scope", "all"),
            "--threads", THREADS,
            "--out-lineages", file.path(outdir, "lineages.tsv"),
            "--out-nodes", file.path(outdir, "nodes.tsv"),
            "--out-node-links", file.path(outdir, "node-links.tsv"),
            "--out-builders", file.path(outdir, "builders.tsv"),
            "--out-aa-cdist", file.path(outdir, "aa-cdist.tsv"),
            "--out-germline-mutations", file.path(outdir, "germline-mutations.tsv"),
            "--out-consensus", file.path(outdir, "consensus.tsv"),
            "--out-anchor-distances", file.path(outdir, "anchor-distances.tsv"))
  log <- system2("Rscript", argv, stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  list(ok = is.null(status), log = paste(c(align_log, log), collapse = "\n"), dir = outdir,
       aligned = read_tsv(aligned),
       lineages = if (is.null(status)) read_tsv(file.path(outdir, "lineages.tsv")),
       nodes = if (is.null(status)) read_tsv(file.path(outdir, "nodes.tsv")),
       links = if (is.null(status)) read_tsv(file.path(outdir, "node-links.tsv")),
       builders = if (is.null(status)) read_tsv(file.path(outdir, "builders.tsv")),
       cdist = if (is.null(status)) read_tsv(file.path(outdir, "aa-cdist.tsv")),
       germline = if (is.null(status)) read_tsv(file.path(outdir, "germline-mutations.tsv")),
       consensus = if (is.null(status)) read_tsv(file.path(outdir, "consensus.tsv")),
       anchors = if (is.null(status)) read_tsv(file.path(outdir, "anchor-distances.tsv")))
}

# All scenarios run in parallel before any check.
RUN_SPECS <- list(
  list("paired", TRUE, "raxml", TRUE),
  list("bulk", TRUE, "raxml", TRUE),
  list("joins", TRUE, "raxml", TRUE),
  list("twins", TRUE, "raxml", TRUE),
  list("table", TRUE, "raxml", FALSE),
  list("marks", TRUE, "raxml", FALSE),
  list("truncated", FALSE, "raxml", FALSE),
  list("gapped", FALSE, "raxml", FALSE),
  list("oddids", TRUE, "raxml", TRUE),
  list("mixedgaps", FALSE, "raxml", FALSE),
  list("tiny", TRUE, "igphyml", TRUE),
  list("empty", TRUE, "raxml", TRUE),
  list("bare", FALSE, "raxml", FALSE))
names(RUN_SPECS) <- vapply(RUN_SPECS, `[[`, character(1), 1)
workers <- max(1L, min(length(RUN_SPECS), parallel::detectCores()))
# One thread per run; the runs already fill the cores.
THREADS <- "1"
started <- Sys.time()
results <- parallel::mclapply(RUN_SPECS, function(spec) {
  run_one(spec[[1]], light = spec[[2]], builder = spec[[3]], airr = spec[[4]])
}, mc.cores = workers, mc.preschedule = FALSE)
cat(sprintf("%d scenarios on %d workers in %.0f s\n", length(RUN_SPECS), workers,
            as.numeric(difftime(Sys.time(), started, units = "secs"))))
run_trees <- function(scenario) {
  value <- results[[scenario]]
  if (inherits(value, "try-error") || is.null(value)) {
    value <- list(ok = FALSE, log = paste("worker failed:", value), dir = NA_character_)
  }
  value
}
tail_of <- function(r) cat(substr(r$log, max(1, nchar(r$log) - 1500), nchar(r$log)), "\n")
# One check over named parts; a failure lists the parts that failed.
ok_all <- function(what, parts) {
  bad <- names(parts)[!vapply(parts, isTRUE, logical(1))]
  ok(if (length(bad)) sprintf("%s [failed: %s]", what, paste(bad, collapse = "; ")) else what,
     length(bad) == 0)
}
# A run that crashed fails its first check and shows the log.
crashed <- function(r, what) {
  ok_all(what, list(runs = FALSE))
  tail_of(r)
}

# aa-cdist comes from the consensus, not the tree, so it is checked against lineage size.
AA_CDIST_MIN_SEQUENCES <- 10
check_aa_cdist <- function(r, label) {
  counts <- setNames(as.integer(r$consensus$consensus_sequence_count),
                     r$consensus$lineage_id)
  scored <- unique(r$lineages$lineage_id[r$lineages$sequence_id %in% r$cdist$sequence_id])
  above <- names(counts)[counts >= AA_CDIST_MIN_SEQUENCES]
  parts <- list(
    # Unscored lineages get a count too.
    "every lineage carries a consensus count" =
      setequal(names(counts), unique(r$lineages$lineage_id[r$lineages$sequence_id %in%
                                                           r$links$sequence_id])) ||
        length(setdiff(names(counts), unique(r$lineages$lineage_id))) == 0,
    "scored exactly at or above the floor" = setequal(scored, above),
    "no sequence scored twice" = !any(duplicated(r$cdist$sequence_id)),
    "scores are non-negative integers" =
      nrow(r$cdist) == 0 || all(!is.na(suppressWarnings(as.integer(r$cdist$aa_cdist))) &
                                as.integer(r$cdist$aa_cdist) >= 0))
  # The consensus is the centre, so some member must sit near it.
  if (nrow(r$cdist)) {
    per <- split(as.integer(r$cdist$aa_cdist),
                 r$lineages$lineage_id[match(r$cdist$sequence_id, r$lineages$sequence_id)])
    parts[["a member near each consensus"]] <- all(vapply(per, function(v) min(v) <= 2, logical(1)))
  }
  ok_all(paste(label, "check_aa_cdist"), parts)
}

# Node sequences: IUPAC, one frame per lineage, heavy on every node, light only where the lineage has light chains.
# With `run`, each chain's width must match its padded, ungapped aligned rows, so the split lost nothing.
IUPAC <- c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-", ".")
check_node_sequences <- function(nodes, label, light = FALSE, run = NULL) {
  what <- paste(label, "check_node_sequences")
  if (!all(c("heavy_sequence", "light_sequence") %in% names(nodes))) {
    return(ok_all(what, list("heavy and light sequence columns" = FALSE)))
  }
  heavy <- nodes$heavy_sequence
  lite <- nodes$light_sequence
  has_light <- !is.na(lite) & nzchar(lite)
  letters_used <- unique(unlist(strsplit(toupper(paste(c(heavy, lite[has_light]), collapse = "")), "")))
  per_lineage <- tapply(nchar(heavy), nodes$lineage_id, function(x) length(unique(x)))
  light_frames <- tapply(ifelse(has_light, nchar(lite), -1L), nodes$lineage_id,
                         function(x) length(unique(x)))
  parts <- list(
    "heavy on every node" = nrow(nodes) == 0 || all(!is.na(heavy) & nchar(heavy) > 0),
    "light only with light chains" = if (light) any(has_light) else !any(has_light),
    "IUPAC only" = all(letters_used %in% IUPAC),
    "one frame per lineage" = length(per_lineage) == 0 || all(per_lineage == 1),
    "heavy is whole codons" = all(nchar(heavy) %% 3L == 0L),
    "light on all nodes of a lineage or none, one frame" =
      length(light_frames) == 0 || all(light_frames == 1))
  if (!is.null(run) && nrow(nodes)) {
    a <- run$aligned
    a$lineage_id <- run$lineages$lineage_id[match(a$sequence_id, run$lineages$sequence_id)]
    a$heavy <- a$locus == "IGH"
    fl <- as.integer(a$frame_left)
    fr <- as.integer(a$frame_right)
    # Pad as trees.R does, then drop columns that are "." in every germline row, as dowser does.
    expected <- function(lid, want_heavy) {
      rows <- which(a$lineage_id %in% lid & a$heavy == want_heavy)
      if (!length(rows)) return(NA_integer_)
      g <- paste0(strrep("N", max(fl[rows]) - fl[rows]), a$germline_alignment[rows],
                  strrep("N", max(fr[rows]) - fr[rows]))
      if (length(unique(nchar(g))) != 1) return(NA_integer_)
      sum(colSums(do.call(rbind, strsplit(g, "")) != ".") > 0)
    }
    lineages <- unique(nodes$lineage_id)
    heavy_width <- vapply(lineages, function(lid) nchar(heavy[nodes$lineage_id == lid][1]), integer(1))
    light_width <- vapply(lineages, function(lid) {
      x <- lite[nodes$lineage_id == lid & has_light]
      if (length(x)) nchar(x[1]) else 0L
    }, integer(1))
    want_heavy <- vapply(lineages, expected, integer(1), want_heavy = TRUE)
    want_light <- vapply(lineages, function(lid) {
      if (!light) return(0L)
      w <- expected(lid, FALSE)
      if (is.na(w) && !any(a$lineage_id %in% lid & !a$heavy)) 0L else w
    }, integer(1))
    parts[["heavy as wide as its rows"]] <-
      all(heavy_width == want_heavy, na.rm = TRUE) && any(!is.na(want_heavy))
    parts[["light as wide as its rows"]] <-
      all(light_width == want_light, na.rm = TRUE) && (!light || any(want_light > 0, na.rm = TRUE))
  }
  ok_all(what, parts)
}

# Depth and per-branch mutations, per chain and alphabet; a chain with no data is blank.
STEP_CHAINS <- c("heavy", "light")
check_node_steps <- function(nodes, label, light = FALSE) {
  what <- paste(label, "check_node_steps")
  for (chain in STEP_CHAINS) for (col in c("mutations_from_parent", "mutation_count_from_parent",
                                           "unresolved_from_parent", "aa_mutations_from_parent",
                                           "aa_mutation_count_from_parent")) {
    column <- paste0(chain, "_", col)
    if (!(column %in% names(nodes))) return(ok_all(what, setNames(list(FALSE), paste("column", column))))
  }
  if (!nrow(nodes)) return(ok_all(what, list()))
  depth <- as.integer(nodes$node_depth)
  parent <- suppressWarnings(as.integer(nodes$parent_id))
  key <- paste(nodes$lineage_id, nodes$node_id)
  depth_of <- setNames(depth, key)
  root <- is.na(parent)
  # Branch fraction of root distance, and tips under the parent (including the node itself).
  frac <- suppressWarnings(as.numeric(nodes$terminal_branch_fraction))
  below <- suppressWarnings(as.integer(nodes$parent_descendant_count))
  observed <- nodes$is_observed == "true"
  parts <- list(
    "one root per lineage at depth zero" =
      all(depth[root] == 0L) && length(unique(nodes$lineage_id[root])) == sum(root),
    "one step below the parent" =
      all(depth[!root] == depth_of[paste(nodes$lineage_id[!root], parent[!root])] + 1L),
    "terminal branch fraction in [0, 1], unset at root" =
      all(is.na(frac[root])) && all(frac[!is.na(frac)] >= 0 & frac[!is.na(frac)] <= 1 + 1e-9),
    "observed node among its parent's descendants" =
      all(is.na(below[root])) && all(below[observed & !root] >= 1L),
    "root acquired nothing" =
      all(nodes$heavy_mutations_from_parent[root] == "") &&
        all(nodes$heavy_aa_mutations_from_parent[root] %in% c("", NA)))
  count_listed <- function(x) ifelse(is.na(x) | x == "", 0L, lengths(strsplit(x, ",")))
  as_int <- function(x) suppressWarnings(as.integer(x))
  for (chain in STEP_CHAINS) {
    nt <- nodes[[paste0(chain, "_mutations_from_parent")]]
    aa <- nodes[[paste0(chain, "_aa_mutations_from_parent")]]
    nt_n <- as_int(nodes[[paste0(chain, "_mutation_count_from_parent")]])
    aa_n <- as_int(nodes[[paste0(chain, "_aa_mutation_count_from_parent")]])
    present <- !is.na(nt_n)
    at <- function(name) paste(chain, name)
    if (chain == "light" && !light) {
      parts[[at("steps absent without light chains")]] <- !any(present)
      next
    }
    parts[[at("steps reported")]] <- if (chain == "light") any(present) else all(present)
    # Entries are base, position, base; only settled bases.
    listed <- unlist(strsplit(nt[present & nt != ""], ","))
    parts[[at("mutations are real changes between settled bases")]] <-
      length(listed) == 0 || (all(grepl("^[ACGT][0-9]+[ACGT]$", listed)) &&
        all(substr(listed, 1, 1) != substr(listed, nchar(listed), nchar(listed))))
    parts[[at("mutation count matches the list")]] <- all(nt_n[present] == count_listed(nt[present]))
    parts[[at("unresolved never negative")]] <-
      all(as_int(nodes[[paste0(chain, "_unresolved_from_parent")]][present]) >= 0)
    # Residue, codon position, residue; never a stop or unreadable codon.
    scored <- present & !is.na(aa_n)
    aa_listed <- unlist(strsplit(aa[scored & aa != ""], ","))
    parts[[at("amino acid changes between settled residues")]] <-
      length(aa_listed) == 0 || all(grepl("^[ACDEFGHIKLMNPQRSTVWY][0-9]+[ACDEFGHIKLMNPQRSTVWY]$", aa_listed))
    parts[[at("amino acid count matches the list")]] <- all(aa_n[scored] == count_listed(aa[scored]))
    parts[[at("residues never exceed bases")]] <- all(aa_n[scored] <= nt_n[scored])
    parts[[at("amino acid figures present")]] <- any(scored) || !any(present)
  }
  ok_all(what, parts)
}

# Every lineage is one rooted tree over contiguous node ids.
check_topology <- function(nodes, label) {
  contiguous <- TRUE; one_root <- TRUE; parents_known <- TRUE
  for (lid in unique(nodes$lineage_id)) {
    d <- nodes[nodes$lineage_id == lid, ]
    ids <- as.integer(d$node_id)
    parents <- suppressWarnings(as.integer(d$parent_id))
    if (!identical(sort(ids), seq_along(ids))) contiguous <- FALSE
    if (sum(is.na(parents) | parents == "") != 1L) one_root <- FALSE
    if (any(!is.na(parents) & !(parents %in% ids))) parents_known <- FALSE
  }
  ok_all(paste(label, "check_topology"),
         list("contiguous node ids" = contiguous, "one root" = one_root,
              "parents exist" = parents_known))
}

# Partial coverage. The frame comes from the junction, so rows must stay whole codons with no germline stop.
# A frame taken from the alignment start would make stops that dowser then filters out.
STOP_CODONS <- c("TAA", "TAG", "TGA")
in_frame_stops <- function(x) {
  s <- strsplit(toupper(x), "")[[1]]
  n <- length(s) %/% 3
  if (!n) return(0L)
  sum(vapply(seq_len(n),
             function(k) paste(s[(3 * k - 2):(3 * k)], collapse = "") %in% STOP_CODONS,
             logical(1)))
}
check_truncated <- function(r, label, junction_only = FALSE) {
  what <- paste(label, "check_truncated")
  if (!isTRUE(r$ok)) return(crashed(r, what))
  a <- r$aligned
  parts <- list(
    "rows survive the rebuild" = nrow(a) > 0,
    "flanks are whole codons" =
      all(as.integer(a$frame_left) %% 3 == 0) && all(as.integer(a$frame_right) %% 3 == 0),
    "rows are whole codons" = all(nchar(a$sequence_alignment) %% 3 == 0),
    # Exporter padding and IMGT gaps both go; a "." left is a deletion against the germline.
    "exporter padding gone" =
      !any(grepl("^\\.", a$sequence_alignment)) && !any(grepl("\\.$", a$sequence_alignment)),
    "no IMGT gap left in the germline" = !any(grepl(".", a$germline_alignment, fixed = TRUE)),
    "no germline in-frame stop" = all(vapply(a$germline_alignment, in_frame_stops, integer(1)) == 0L),
    "tree step keeps its sequences" = !grepl("No clones remain after makeAirrClone", r$log))
  if (junction_only) {
    parts[["no flank left"]] <- all(as.integer(a$frame_left) == 0L) && all(as.integer(a$frame_right) == 0L)
    parts[["germline all N"]] <- all(grepl("^N+$", a$germline_alignment))
  } else {
    parts[["trees are built"]] <- nrow(r$nodes) > 0
  }
  ok_all(what, parts)
  if (!junction_only) check_topology(r$nodes, label)
}

cat("== paired: light chains, split, dark, anchor, insertions, a slash in the ids ==\n")
p <- run_trees("paired")
first <- "paired: every clonotype gets one lineage, both chains joined"
if (p$ok) {
  here <- file.path(root, "paired")
  clones <- read_tsv(file.path(here, "clones.tsv"))
  lightless <- readLines(file.path(here, "light-less.txt"))
  ok_all(first, list(
    "one lineage each" = setequal(p$lineages$sequence_id, clones$sequence_id) &&
      !any(duplicated(p$lineages$sequence_id)),
    "no blank id" = all(nzchar(p$lineages$lineage_id)) && !any(is.na(p$lineages$lineage_id)),
    "light chain joined" = grepl("matched \\d+ clonotypes to a light chain", p$log)))
  # Key property: a light-less clonotype is not held apart from its paired kin.
  ll <- p$lineages[p$lineages$sequence_id %in% lightless, ]
  paired_lineages <- p$lineages$lineage_id[!(p$lineages$sequence_id %in% lightless)]
  ok("paired: light-less clonotypes share lineages with paired ones",
     nrow(ll) > 0 && any(ll$lineage_id %in% paired_lineages))

  target <- readLines(file.path(here, "split-clone.txt"))
  recombined <- readLines(file.path(here, "recombined.txt"))
  members <- p$lineages[startsWith(p$lineages$lineage_id, paste0(target, "_")), ]
  ok_all("paired: the light V/J split separates the recombined members", list(
    "clone splits" = length(unique(members$lineage_id)) > 1,
    "recombined apart" = length(intersect(members$lineage_id[members$sequence_id %in% recombined],
                                          members$lineage_id[!(members$sequence_id %in% recombined)])) == 0))

  # RAxML fails to partition a one-chain clone, and getTrees would silently drop it.
  dark <- readLines(file.path(here, "dark.txt"))
  dark_lineages <- unique(p$lineages$lineage_id[p$lineages$sequence_id %in% dark])
  built <- unique(p$nodes$lineage_id)
  ok_all("paired: lineages with no light chain still get trees", list(
    "dark built" = length(dark_lineages) == 2 && all(dark_lineages %in% built),
    "paired built" = length(setdiff(built, dark_lineages)) > 0))
  # Two tips plus germline is one short of RAxML-NG's minimum.
  builder_of <- setNames(p$builders$tree_builder, p$builders$lineage_id)
  tips <- table(p$lineages$lineage_id)
  ok_all("paired: two tips by parsimony, larger lineages by FastTree+RAxML", list(
    "pratchet" = any(tips == 2) && all(builder_of[names(tips)[tips == 2]] == "pratchet"),
    "fasttree-raxml" = all(builder_of[names(tips)[tips >= 3]] == "fasttree-raxml")))
  ok_all("paired: tree build partitions per chain and collapses internal nodes", list(
    "RAxML partitions" = grepl("partition: scaled", p$log),
    "collapse ran" = grepl("collapsed \\d+ of \\d+ internal nodes", p$log) &&
      !grepl("node collapse failed", p$log)))

  # dowser puts the lineage id in a temp path and an unquoted shell word.
  ok_all("paired: a slash in ids loses no clone and ids come back", list(
    "no clone lost" = !grepl("Tree building failed", p$log),
    "ids come back" = all(grepl("^mouse D0_1/", unique(p$nodes$lineage_id))) &&
      all(unique(p$nodes$lineage_id) %in% p$lineages$lineage_id)))

  inserted <- as.integer(readLines(file.path(here, "inserted-rows.txt")))
  ok_all("paired: insertions dropped, no lineage skipped for length or stops", list(
    "dropped and counted" = grepl(sprintf("insertion columns dropped from %d of \\d+ rows", inserted), p$log),
    "no length skip" = !grepl("alignment lengths differ", p$log),
    "no stop removal" = !grepl("inframe stop codon", p$log)))

  check_topology(p$nodes, "paired")
  check_node_sequences(p$nodes, "paired", light = TRUE, run = p)
  check_node_steps(p$nodes, "paired", light = TRUE)
  ok_all("paired: node links are valid and every observed node is linked", list(
    "links match lineages" = all(paste(p$links$sequence_id, p$links$lineage_id) %in%
                                   paste(p$lineages$sequence_id, p$lineages$lineage_id)),
    "observed nodes linked" = sum(p$nodes$is_observed == "true") == nrow(p$links)))

  # A light figure needs a light chain at both ends.
  anchor_id <- readLines(file.path(here, "anchor.txt"))
  d <- p$anchors
  as_int <- function(x) suppressWarnings(as.integer(x))
  no_light <- d$sequence_id %in% lightless
  with_light <- !no_light & !is.na(as_int(d$anchor_nt_light))
  ok_all("paired: anchor distances, heavy set, light only with both chains", list(
    "relatives to the anchor" = nrow(d) > 0 && all(d$anchor_id == anchor_id) && !(anchor_id %in% d$sequence_id),
    "heavy set, aa within nt" = all(!is.na(as_int(d$anchor_nt_heavy))) &&
      all(as_int(d$anchor_aa_heavy) <= as_int(d$anchor_nt_heavy)),
    "no light figure without light" = any(no_light) && all(is.na(as_int(d$anchor_nt_light[no_light]))),
    "light figure when paired, aa within nt" = any(with_light) &&
      all(as_int(d$anchor_aa_light[with_light]) <= as_int(d$anchor_nt_light[with_light]))))
} else crashed(p, first)

cat("== twins: one heavy chain with two light chains votes once ==\n")
tw <- run_trees("twins")
first <- "twins: a twin shares its source's lineage and score and adds no vote"
if (!isTRUE(tw$ok)) crashed(tw, first) else {
  here <- file.path(root, "twins")
  twins <- readLines(file.path(here, "twins.txt"))
  sources <- readLines(file.path(here, "sources.txt"))
  lineage_of <- setNames(tw$lineages$lineage_id, tw$lineages$sequence_id)
  scored <- setNames(tw$cdist$aa_cdist, tw$cdist$sequence_id)
  counts <- setNames(as.integer(tw$consensus$consensus_sequence_count), tw$consensus$lineage_id)
  lid <- unique(lineage_of[sources])
  # Every other member of the deep clone carries a distinct heavy sequence.
  members <- names(lineage_of)[lineage_of %in% lid]
  ok_all(first, list(
    "same lineage" = length(lid) == 1 && identical(unname(lineage_of[twins]), unname(lineage_of[sources])),
    "scored" = all(c(twins, sources) %in% names(scored)),
    "same score" = identical(unname(scored[twins]), unname(scored[sources])),
    "one vote per heavy sequence" = length(lid) == 1 &&
      counts[[lid]] == length(members) - sum(members %in% twins)))
  check_aa_cdist(tw, "twins")
}

cat("== bulk: 5'-truncated light-less members join the light subgroup of their cell ==\n")
b <- run_trees("bulk")
first <- "bulk: each truncated bulk copy lands in its source cell's lineage"
if (b$ok) {
  here <- file.path(root, "bulk")
  bulk <- readLines(file.path(here, "bulk.txt"))
  lineage_of <- setNames(b$lineages$lineage_id, b$lineages$sequence_id)
  source_of <- sub("^ck_bulk_", "ck_", bulk)
  ok_all(first, list(
    "clone splits" = length(unique(lineage_of[source_of])) > 1,
    "with its source" = all(lineage_of[bulk] == lineage_of[source_of])))
} else crashed(b, first)

cat("== joins: shared V/J/junction, one light chain across many clonotypes ==\n")
j <- run_trees("joins")
first <- "joins: every clonotype reaches its heavy and light chain and gets a lineage"
if (j$ok) {
  n <- nrow(read_tsv(file.path(root, "joins", "clonotypes.tsv")))
  # Joining the wrong way round would give a row to one clonotype and drop the rest.
  ok_all(first, list(
    "heavy" = grepl(sprintf("matched %d clonotypes to a heavy chain", n), j$log),
    "light" = grepl(sprintf("matched %d clonotypes to a light chain", n), j$log),
    "lineage" = length(unique(j$lineages$sequence_id)) == n))
} else crashed(j, first)

cat("== table: heavy only, alignments on the table, aa-cdist, copies, anchor ==\n")
tb <- run_trees("table")
first <- "table: light chain resolution off, lineages are the heavy clones"
if (tb$ok) {
  here <- file.path(root, "table")
  clones <- read_tsv(file.path(here, "clones.tsv"))
  ok_all(first, list(
    "resolution off" = grepl("light chain resolution off", tb$log),
    "heavy clones" = all(tb$lineages$lineage_id %in% clones$clone_id)))
  # Rows are keyed by clonotype id, so the ambiguous tuple join must not run.
  from_table <- regmatches(tb$log, regexec("clonotype table: (\\d+) heavy chains", tb$log))[[1]]
  matched <- regmatches(tb$log, regexec("matched (\\d+) clonotypes to a heavy chain", tb$log))[[1]]
  ok("table: alignments come off the table and the AIRR join is not attempted",
     length(from_table) == 2 && length(matched) == 2 && from_table[2] == matched[2] &&
       as.integer(matched[2]) == nrow(clones))
  ok("table: RAxML node sequences are thresholded to a credible set",
     !grepl("ancestral sequences kept unthresholded", tb$log) &&
       !grepl("keeping RAxML's own node sequences", tb$log))
  check_topology(tb$nodes, "table")
  check_node_sequences(tb$nodes, "table", run = tb)
  check_node_steps(tb$nodes, "table")

  counts <- setNames(as.integer(tb$consensus$consensus_sequence_count), tb$consensus$lineage_id)
  v <- as.integer(tb$cdist$aa_cdist)
  ok_all("table: aa-cdist scores a lineage above the floor and discriminates", list(
    "above the floor" = any(counts >= AA_CDIST_MIN_SEQUENCES) && nrow(tb$cdist) > 0,
    "scores differ" = length(unique(v)) > 1))
  check_aa_cdist(tb, "table")
  g <- suppressWarnings(as.integer(tb$germline$germline_mutation_count))
  ok("table: germline mutations, one non-negative count per aligned clonotype, some mutated",
     nrow(tb$germline) > 0 && !any(duplicated(tb$germline$sequence_id)) &&
       all(!is.na(g) & g >= 0) && any(g > 0) &&
       all(tb$germline$sequence_id %in% tb$aligned$sequence_id))

  duplicated_ids <- readLines(file.path(here, "duplicated.txt"))
  anchor_id <- readLines(file.path(here, "anchor.txt"))
  groups <- setNames(tb$lineages$group_id, tb$lineages$sequence_id)
  nodes_of <- tb$links[tb$links$sequence_id %in% duplicated_ids, ]
  rep_of <- nodes_of$sequence_id[nodes_of$is_representative == "true"]
  distinct <- tapply(tb$lineages$group_id, tb$lineages$lineage_id, function(g) length(unique(g)))
  scored <- setNames(tb$cdist$aa_cdist, tb$cdist$sequence_id)
  ok_all("table: copies group, link, represent, vote once and share a score", list(
    "grouping" = length(unique(groups[duplicated_ids])) == 2 && sum(table(tb$lineages$group_id) > 1) == 2,
    "links" = setequal(nodes_of$sequence_id, duplicated_ids) &&
      length(unique(paste(nodes_of$lineage_id, nodes_of$node_id))) == 2,
    "representative" = length(rep_of) == 2 && anchor_id %in% rep_of,
    "one vote per sequence" = all(counts[names(distinct)] == distinct),
    "shared score" = all(tapply(scored[names(groups)[names(groups) %in% names(scored)]],
                                groups[names(groups) %in% names(scored)],
                                function(x) length(unique(x))) == 1)))
  ok("table: heavy-only anchor relatives have no light figure",
     nrow(tb$anchors) > 0 && all(tb$anchors$anchor_id == anchor_id) &&
       all(is.na(suppressWarnings(as.integer(tb$anchors$anchor_nt_light)))))
} else crashed(tb, first)

cat("== truncated: full, FR2, CDR2 and '.'-padded coverage in one lineage ==\n")
tr <- run_trees("truncated")
check_truncated(tr, "truncated")
if (isTRUE(tr$ok)) {
  ok_all("truncated: coverages differ and short members are padded", list(
    "coverages differ" = length(unique(as.integer(tr$aligned$frame_left))) > 2,
    "padded" = grepl("padded to a common frame", tr$log) && !grepl("alignment lengths differ", tr$log)))
}

cat("== marks: anchors read from merge's annotations table ==\n")
mk <- run_trees("marks")
tb2 <- run_trees("table")
if (isTRUE(mk$ok) && isTRUE(tb2$ok)) {
  key <- function(a) sort(paste(a$sequence_id, a$anchor_id))
  ok("marks: anchors come from the annotations table, as they did from the clonotype table",
     nrow(mk$anchors) > 0 && identical(key(mk$anchors), key(tb2$anchors)))
} else crashed(if (isTRUE(mk$ok)) tb2 else mk, "marks: runs")

cat("== oddids: ':', ';', ',', '=' and spaces in clonotype ids ==\n")
od <- run_trees("oddids")
if (isTRUE(od$ok)) {
  clones <- read_tsv(file.path(root, "oddids", "clones.tsv"))
  ok_all("oddids: trees are built and every tip links to its own clonotype id", list(
    "trees are built" = nrow(od$nodes) > 0,
    "links resolve" = nrow(od$links) > 0 && all(od$links$sequence_id %in% clones$sequence_id),
    "ids kept as given" = any(grepl("k:a;b,c=d e", od$nodes$label, fixed = TRUE))))
} else crashed(od, "oddids: runs")

cat("== mixedgaps: IMGT-gapped and ungapped rows in one donor ==\n")
gp <- run_trees("gapped")
mg <- run_trees("mixedgaps")
if (isTRUE(gp$ok) && isTRUE(mg$ok)) {
  by_id <- function(a) a[order(a$sequence_id), c("sequence_id", "sequence_alignment", "germline_alignment")]
  ok_all("mixedgaps: gaps removed, so both conventions rebuild to the same rows and trees", list(
    "says so" = grepl("IMGT gaps removed from", mg$log, fixed = TRUE),
    "same rows as all gapped" = identical(by_id(gp$aligned), by_id(mg$aligned)),
    "no gap left in a germline" = !any(grepl(".", mg$aligned$germline_alignment, fixed = TRUE)),
    "trees are built" = nrow(mg$nodes) > 0 && nrow(mg$nodes) == nrow(gp$nodes)))
} else crashed(if (isTRUE(gp$ok)) mg else gp, "mixedgaps: runs")

cat("== tiny: IgPhyML, two tips per lineage ==\n")
g <- run_trees("tiny")
first <- "tiny: IgPhyML builds the trees and is recorded as the builder"
if (g$ok) {
  ok(first, nrow(g$nodes) > 0 && nrow(g$builders) > 0 && all(g$builders$tree_builder == "igphyml"))
  ok("tiny: IgPhyML gets one omega per chain", grepl("partition: hl", g$log))
  check_topology(g$nodes, "tiny")
  check_node_sequences(g$nodes, "tiny", light = TRUE, run = g)
  check_node_steps(g$nodes, "tiny", light = TRUE)
} else crashed(g, first)

cat("== empty: a donor that contributed nothing ==\n")
e <- run_trees("empty")
first <- "empty: runs, writes empty tables and says why"
if (e$ok) {
  ok_all(first, list(
    "empty tables" = nrow(e$lineages) == 0 && nrow(e$nodes) == 0 && nrow(e$links) == 0,
    # The UI parses this exact text.
    "says why" = grepl("no trees for this group", e$log)))
} else crashed(e, first)

cat("== bare: no alignments anywhere ==\n")
ba <- run_trees("bare")
first <- "bare: runs, every clonotype keeps its lineage with no tree, the reason names the missing input"
if (ba$ok) {
  clones <- read_tsv(file.path(root, "bare", "clones.tsv"))
  ok_all(first, list(
    "lineages kept, no tree" = setequal(ba$lineages$sequence_id, clones$sequence_id) && nrow(ba$nodes) == 0,
    "reason" = grepl("no clonotype carries a heavy chain alignment: no trees for this group", ba$log)))
} else crashed(ba, first)

cat(sprintf("\n%d checks, %d failures\n", checks, failures))
quit(status = if (failures > 0L) 1L else 0L)
