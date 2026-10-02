#!/usr/bin/env Rscript
# Two stages, chosen with `--stage`. align: joins each clonotype to its alignment and
# rebuilds it in its germline frame, before clustering, so HILARy reads the same mutations.
# trees: splits HILARy's heavy-only clones by light chain into final lineages, builds trees.

# Libraries load in the trees stage only; attaching dowser costs ~19 s that align does not need.

# Print each warning as it happens: dowser drops a failed clone with only a warning.
options(warn = 1)

# Send warnings to stdout, since the log sink does not capture the message stream.
globalCallingHandlers(warning = function(w) {
  cat(sprintf("warning: %s\n", conditionMessage(w)))
  invokeRestart("muffleWarning")
})

local({
  me <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(normalizePath(me)), "common.R"))
})
stage <- if ("--stage" %in% args) args[match("--stage", args) + 1] else "trees"

# MiXCR exports, `ds<dataset>__<sample>.tsv`; absent when no dataset has a clns.
airr_dir <- opt("--airr-dir", required = FALSE)
# The align stage's output and the trees stage's input.
aligned_path <- opt("--out-aligned", required = FALSE)
aligned_in_path <- opt("--aligned", required = FALSE)
clones_path <- opt("--clones", required = stage == "trees")
clonotypes_path <- opt("--clonotypes")
lineages_path <- opt("--out-lineages", required = stage == "trees")
nodes_path <- opt("--out-nodes", required = stage == "trees")
links_path <- opt("--out-node-links", required = stage == "trees")
builders_path <- opt("--out-builders", required = FALSE)
cdist_path <- opt("--out-aa-cdist", required = FALSE)
germline_mutations_path <- opt("--out-germline-mutations", required = FALSE)
support_path <- opt("--out-consensus", required = FALSE)
# Our wrapper, not raxml-ng: FastTree drafts the topology, as ML search fails past ~1,000 tips.
# The image sets these; natively the programs are on PATH.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
on_path <- function(name, fallback) { p <- unname(Sys.which(name)); if (nzchar(p)) p else fallback }
raxml <- Sys.getenv("RAXML_EXEC", file.path(script_dir, "raxml-ng-fasttree"))
# Windows has no bash or fork: run the wrapper's R port through a batch file. The name
# keeps the "-ng" dowser checks for; short paths, since dowser does not quote the command.
ON_WINDOWS <- .Platform$OS.type == "windows"
if (ON_WINDOWS && !nzchar(Sys.getenv("RAXML_EXEC"))) {
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  raxml <- file.path(tempdir(), "raxml-ng-fasttree.cmd")
  writeLines(c(sprintf('@"%s" "%s" %%*', utils::shortPathName(rscript),
                       utils::shortPathName(file.path(script_dir, "raxml-ng-fasttree.R"))),
               "@exit /b %ERRORLEVEL%"), raxml)
  raxml <- utils::shortPathName(raxml)
}
# Dowser's own nproc forks (mclapply), which Windows cannot; one process there.
fork_workers <- function(n) if (ON_WINDOWS) 1L else n
igphyml <- Sys.getenv("IGPHYML_EXEC", on_path("igphyml", "/usr/local/share/igphyml/src/igphyml"))
# Native igphyml finds its hotspot tables only through IGPHYML_PATH; the image has them at its built-in path.
motifs <- file.path(dirname(igphyml), "..", "share", "igphyml", "motifs")
if (is.na(Sys.getenv("IGPHYML_PATH", NA)) && dir.exists(motifs)) Sys.setenv(IGPHYML_PATH = normalizePath(motifs))
# Which lineages IgPhyML builds: "none", "anchored" or "all". FastTree+RAxML build the rest.
igphyml_scope <- if (is.null(opt("--igphyml-scope", required = FALSE))) "none" else opt("--igphyml-scope")
anchor_path <- opt("--out-anchor-distances", required = FALSE)
# A lineage above this many tips is subsampled to it. Absent means no cap.
max_tips <- {
  m <- opt("--max-tips", required = FALSE)
  if (is.null(m)) Inf else max(3L, as.integer(m))
}
# A lineage below this many tips keeps its membership and gets no tree.
min_tips <- {
  m <- opt("--min-tips", required = FALSE)
  if (is.null(m)) 2L else max(2L, as.integer(m))
}
# Off for bulk or when the user declines it; turned off below if no light rows arrive.
use_light <- has_flag("--light-chains")
# RAxML runs one process per lineage; IgPhyML uses this as its thread count.
threads <- {
  t <- opt("--threads", required = FALSE)
  if (is.null(t)) 1L else max(1L, as.integer(t))
}

# Elapsed wall time as H:MM:SS, from the start of this run unless told otherwise.
STARTED <- Sys.time()
clock <- function(since = STARTED) {
  s <- as.integer(round(as.numeric(difftime(Sys.time(), since, units = "secs"))))
  sprintf("%d:%02d:%02d", s %/% 3600L, s %/% 60L %% 60L, s %% 60L)
}

# Tee output to a per-donor log file; `split = TRUE` keeps stdout for the platform log.
log_path <- opt("--out-log", required = FALSE)
if (!is.null(log_path)) sink(file(log_path, open = "wt"), split = TRUE)
close_log <- function() {
  cat(sprintf("%s stage finished in %s\n", stage, clock()))
  while (sink.number() > 0) sink()
}

# Forked workers that fit in memory, for dowser's own forks. GC makes a fork drift to a
# full copy of the parent, so allow 80% of the cgroup limit over the parent's size.
memory_limit <- function() {
  for (path in c("/sys/fs/cgroup/memory.max",
                 "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    if (!file.exists(path)) next
    v <- suppressWarnings(as.numeric(readLines(path, n = 1, warn = FALSE)))
    # An unset cgroup limit reads as "max" or as a number near the word size.
    if (!is.na(v) && v > 0 && v < 2^60) return(v)
  }
  NA_real_
}
self_rss <- function() {
  if (!file.exists("/proc/self/statm")) return(NA_real_)
  v <- suppressWarnings(as.numeric(strsplit(readLines("/proc/self/statm", n = 1, warn = FALSE), " ")[[1]][2]))
  if (is.na(v)) NA_real_ else v * 4096
}
fit_workers <- function(want, what) {
  limit <- memory_limit()
  rss <- self_rss()
  if (is.na(limit) || is.na(rss) || rss <= 0) {
    cat(sprintf("%s: %d workers (no memory limit readable, so not bounded)\n", what, want))
    return(want)
  }
  fits <- max(1L, as.integer(floor(limit * 0.8 / rss)) - 1L)
  got <- min(want, fits)
  cat(sprintf("%s: %d of %d workers (%.1f GiB limit, %.1f GiB resident before forking)\n",
              what, got, want, limit / 2^30, rss / 2^30))
  got
}


HEAVY <- "IGH"

# Progress prefix; the block reads the last line into a bar (tips done over all tips).
PROGRESS_PREFIX <- "[==PROGRESS==]"
# Stamped with the elapsed time up front, since the UI reads the percentage off the end.
progress <- function(text) cat(sprintf("%s [%s] %s\n", PROGRESS_PREFIX, clock(), text))
pct <- function(n, total) if (total > 0) 100 * n / total else 0
# The bar moves after every finished lineage, at most this often, in seconds.
TREE_PROGRESS_EVERY <- 1
trees_total <- 0
trees_done <- 0
last_tree_progress <- Sys.time()

# exportAirr writes no locus column, so derive it from the V gene.
locus_of <- function(v_call) substr(gene_of(v_call), 1, 3)

# Input.

# A run is one donor, so an empty result writes empty tables (finish_empty), not an error.
ALIGNED_COLUMNS <- c("sequence_id", "v_call", "j_call", "junction", "sequence_alignment",
                     "germline_alignment", "masked_sequence_alignment", "locus",
                     "frame_left", "frame_right")

# Per-branch changes into a node, by chain and alphabet; light columns blank without light.
STEP_COLUMNS <- unlist(lapply(c("heavy", "light"), function(chain) paste0(chain, c(
  "_mutations_from_parent", "_mutation_count_from_parent", "_unresolved_from_parent",
  "_aa_mutations_from_parent", "_aa_mutation_count_from_parent"))))
NODE_COLUMNS <- c("lineage_id", "node_id", "parent_id", "distance", "is_observed", "label",
                  "heavy_sequence", "light_sequence", "node_depth",
                  "terminal_branch_fraction", "parent_descendant_count", STEP_COLUMNS)
# Every distance is to the one anchor `anchor_id` names, chosen by the heavy chain.
# Heavy chain differences from the germline, per clonotype.
GERMLINE_MUTATION_COLUMNS <- c("sequence_id", "lineage_id", "germline_mutation_count")
ANCHOR_COLUMNS <- c("sequence_id", "anchor_id", "anchor_aa_heavy", "anchor_nt_heavy",
                    "anchor_aa_light", "anchor_nt_light")

finish_empty <- function(reason, membership = NULL) {
  if (stage == "align") {
    # No alignments is a normal state; the trees stage reports it.
    cat(sprintf("%s: no alignments for this group\n", reason))
    empty_aligned <- as.data.frame(setNames(replicate(length(ALIGNED_COLUMNS), character(0), simplify = FALSE), ALIGNED_COLUMNS))
    write.table(empty_aligned, aligned_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    close_log()
    quit(status = 0)
  }
  cat(sprintf("%s: no trees for this group\n", reason))
  if (is.null(membership)) {
    clones_only <- read_tsv(clones_path)
    # No alignments reached these, so each clonotype is its own sequence group.
    membership <- data.frame(sequence_id = clones_only$sequence_id,
                             lineage_id = clones_only$clone_id,
                             group_id = clones_only$sequence_id,
                             link = 1L, stringsAsFactors = FALSE)
  }
  write.table(membership, lineages_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  empty <- function(cols) {
    d <- as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols))
    d
  }
  write.table(empty(NODE_COLUMNS), nodes_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(empty(c("lineage_id", "node_id", "sequence_id", "is_representative", "link")),
              links_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!is.null(builders_path)) {
    write.table(empty(c("lineage_id", "tree_builder")),
                builders_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  }
  # aa-cdist needs no tree, so it may already be written by the time we exit here.
  write_absent <- function(path, cols) {
    if (!is.null(path) && !file.exists(path)) {
      write.table(empty(cols), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    }
  }
  write_absent(cdist_path, c("sequence_id", "lineage_id", "aa_cdist"))
  write_absent(germline_mutations_path, GERMLINE_MUTATION_COLUMNS)
  write_absent(support_path, c("lineage_id", "consensus_sequence_count"))
  write_absent(anchor_path, ANCHOR_COLUMNS)
  close_log()
  quit(status = 0)
}

if (stage == "trees") {
  clones_preview <- read_tsv(clones_path)
  if (!nrow(clones_preview)) finish_empty("no clonotypes in this group",
                                          data.frame(sequence_id = character(0),
                                                     lineage_id = character(0),
                                                     group_id = character(0),
                                                     link = integer(0)))
}

clono <- read_tsv(clonotypes_path)
clones <- if (stage == "trees") read_tsv(clones_path) else NULL

has_light_columns <- all(c("v_call_light", "j_call_light", "junction_light") %in% names(clono))
if (use_light && !has_light_columns) {
  cat("no light chain columns in the clonotype table: light chain resolution off\n")
  use_light <- FALSE
}

present <- function(x) !is.na(x) & nzchar(x)

# Anchors, as `merge` marked them; the "anchored" IgPhyML scope builds their lineages.
anchor_ids <- if ("is_anchor" %in% names(clono)) clono$sequence_id[clono$is_anchor == "true"] else character(0)

if (stage == "align") {
progress("Joining alignments to clonotypes")
# Two alignment sources, maybe both in one run: MiXCR exports under `--airr-dir`, joined
# by V, J and junction, and alignment columns an imported dataset carries, keyed by id.

ALIGNMENT_COLUMNS <- c("sequence_alignment", "germline_alignment")
AIRR_REQUIRED <- c("v_call", "j_call", "junction", ALIGNMENT_COLUMNS)
CHAIN_COLUMNS <- c("sequence_id", "v_call", "j_call", "junction", ALIGNMENT_COLUMNS, "locus")

read_airr <- function(files) {
  if (!length(files)) return(NULL)
  # Only the columns the join and the frame read; exportAirr writes some thirty more.
  keep <- c("sequence_id", AIRR_REQUIRED)
  airr <- do.call(rbind, lapply(files, function(f) {
    airr <- read_tsv(f)
    for (col in AIRR_REQUIRED) {
      if (!col %in% names(airr)) stop("AIRR table is missing ", col)
    }
    airr[, intersect(keep, names(airr)), drop = FALSE]
  }))
  airr$locus <- locus_of(airr$v_call)
  airr$v_gene <- gene_of(airr$v_call)
  airr$j_gene <- gene_of(airr$j_call)
  airr
}

# One entry per dataset, from the `ds<dataset>__` prefix; no prefix means all clonotypes.
tuple_airr <- list()
if (!is.null(airr_dir)) {
  files <- list.files(airr_dir, pattern = "\\.tsv$", full.names = TRUE)
  names_ <- basename(files)
  dataset_of_file <- ifelse(grepl("^ds[^_]+__", names_), sub("^ds([^_]+)__.*$", "\\1", names_), "")
  for (ds in unique(dataset_of_file)) {
    tuple_airr[[ds]] <- read_airr(files[dataset_of_file == ds])
  }
}

clono_dataset <- if ("dataset" %in% names(clono)) clono$dataset else rep("", nrow(clono))
in_dataset <- function(ds) if (nzchar(ds)) clono_dataset == ds else rep(TRUE, nrow(clono))

# AIRR rows to clonotype ids, per locus. Row ids are per file and per cell, so join on
# V gene, J gene and junction within one dataset; the assembled sequence breaks ties.
# The AIRR cell id does not match the block's key, so the clonotype key is the cell id.

assign_ids <- function(rows, keys, main_seqs, ids, label) {
  # rows: AIRR rows of one locus. keys/main_seqs/ids: the clonotype side.
  if (!nrow(rows) || !length(ids)) return(rows[0, , drop = FALSE])
  # A cell-tagged clns repeats each clone per cell; collapse copies, sorted so the kept one is stable.
  before <- nrow(rows)
  rows <- rows[order(rows$v_gene, rows$j_gene, rows$junction, rows$sequence_alignment), ]
  rows <- rows[!duplicated(rows[, c("v_gene", "j_gene", "junction", "sequence_alignment")]), ]
  if (nrow(rows) < before) {
    cat(sprintf("%s: %d of %d AIRR rows were per-cell copies of another\n",
                label, before - nrow(rows), before))
  }

  ungapped <- toupper(gsub("[.-]", "", rows$sequence_alignment))
  row_keys <- paste(rows$v_gene, rows$j_gene, rows$junction, sep = "\r")

  # Per clonotype, not per row: one chain's clone can pair into many clonotypes.
  idx <- split(seq_along(row_keys), row_keys)
  bucket <- match(keys, names(idx))
  pick <- rep(NA_integer_, length(ids))
  for (i in seq_along(ids)) {
    if (is.na(bucket[i])) next
    cand <- idx[[bucket[i]]]
    if (length(cand) > 1 && !is.na(main_seqs[i]) && nzchar(main_seqs[i])) {
      # Tie: clones assembled wider than CDR3; the assembled sequence separates them.
      hit <- cand[grepl(toupper(main_seqs[i]), ungapped[cand], fixed = TRUE)]
      if (length(hit)) cand <- hit
    }
    pick[i] <- cand[1]
  }

  matched <- !is.na(pick)
  if (any(!matched)) {
    cat(sprintf("%s: %d of %d clonotypes had no AIRR row (%.1f%%)\n",
                label, sum(!matched), length(ids), 100 * sum(!matched) / length(ids)))
  }
  out <- rows[pick[matched], , drop = FALSE]
  out$sequence_id <- ids[matched]
  out
}

# Clonotype table alignments, for clonotypes no AIRR row reached; missing either means no tree.
from_table <- function(suffix, eligible, taken, label) {
  cols <- paste0(ALIGNMENT_COLUMNS, suffix)
  if (!all(cols %in% names(clono))) return(NULL)
  ok <- eligible & present(clono[[cols[1]]]) & present(clono[[cols[2]]]) &
    !(clono$sequence_id %in% taken)
  if (!any(ok)) return(NULL)
  v <- clono[[paste0("v_call", suffix)]][ok]
  out <- data.frame(
    sequence_id = clono$sequence_id[ok],
    v_call = v, j_call = clono[[paste0("j_call", suffix)]][ok],
    junction = clono[[paste0("junction", suffix)]][ok],
    sequence_alignment = clono[[cols[1]]][ok],
    germline_alignment = clono[[cols[2]]][ok],
    locus = if (nzchar(suffix)) locus_of(v) else HEAVY, stringsAsFactors = FALSE)
  cat(sprintf("alignments read from the clonotype table: %d %s chains\n", nrow(out), label))
  out
}

gather <- function(is_locus, eligible, keys, main_seqs, suffix, label) {
  parts <- list()
  # By position: the unprefixed dataset is named "", and `[[""]]` matches nothing.
  for (k in seq_along(tuple_airr)) {
    ds <- names(tuple_airr)[k]
    airr <- tuple_airr[[k]]
    if (is.null(airr)) next
    mine <- eligible & in_dataset(ds)
    rows <- airr[is_locus(airr$locus), , drop = FALSE]
    tag <- if (nzchar(ds)) paste0("dataset ", ds, " ", label) else label
    parts[[length(parts) + 1]] <- assign_ids(rows, keys[mine], main_seqs[mine],
                                             clono$sequence_id[mine], tag)
  }
  taken <- unlist(lapply(parts, function(p) p$sequence_id))
  parts[[length(parts) + 1]] <- from_table(suffix, eligible, taken, label)
  parts <- Filter(function(p) !is.null(p) && nrow(p) > 0, parts)
  if (!length(parts)) return(NULL)
  do.call(rbind, lapply(parts, function(p) p[, CHAIN_COLUMNS, drop = FALSE]))
}

heavy_keys <- paste(gene_of(clono$v_call), gene_of(clono$j_call), clono$junction, sep = "\r")
heavy_main <- if ("main_sequence" %in% names(clono)) clono$main_sequence else rep(NA_character_, nrow(clono))
joined <- gather(function(l) l == HEAVY, rep(TRUE, nrow(clono)), heavy_keys, heavy_main, "", "heavy")
if (is.null(joined)) finish_empty("no clonotype carries a heavy chain alignment")
cat(sprintf("matched %d clonotypes to a heavy chain\n", nrow(joined)))

if (use_light) {
  light_present <- present(clono$v_call_light)
  light_keys <- paste(gene_of(clono$v_call_light), gene_of(clono$j_call_light),
                      clono$junction_light, sep = "\r")
  light_main <- if ("main_sequence_light" %in% names(clono)) clono$main_sequence_light else rep(NA_character_, nrow(clono))
  light <- if (any(light_present)) {
    gather(function(l) l != HEAVY, light_present, light_keys, light_main, "_light", "light")
  } else {
    NULL
  }
  if (is.null(light)) {
    cat("no light chain alignments matched: light chain resolution off\n")
    use_light <- FALSE
  } else {
    cat(sprintf("matched %d clonotypes to a light chain\n", nrow(light)))
    joined <- rbind(joined, light)
  }
}

# The AIRR tables and the clonotype table are done with; the frame below is the peak.
rm(tuple_airr, clono)
invisible(gc())

# Rebuild every row in its germline's frame so a lineage's members share one length:
# V side to the junction, J side after it, query insertions dropped (logged), and the
# junction from the query over an all-N germline (this also masks D, whose stops igphyml
# refuses). Each flank stops at the last real base, not "." padding, and is cut to whole
# codons from the junction. Whole-codon insertions are first slid to a codon boundary,
# since cut mid-codon they make false codons, often stops.

# Works on raw bytes rather than single characters: paste(collapse) was most of the cost.
DOT <- charToRaw("."); DASH <- charToRaw("-"); NB <- charToRaw("N")

# Slides each whole-codon insertion run in `cols` of `g` back to a codon boundary.
codon_align_insertions <- function(g, cols) {
  if (!length(cols)) return(g)
  ins <- cols[g[cols] == DASH]
  if (!length(ins)) return(g)
  breaks <- which(diff(ins) != 1)
  starts <- ins[c(1, breaks + 1)]
  ends <- ins[c(breaks, length(ins))]
  for (r in seq_along(starts)) {
    len <- ends[r] - starts[r] + 1
    if (len %% 3 != 0) next
    before <- cols[cols < starts[r]]
    k <- sum(g[before] != DASH) %% 3
    if (k == 0 || starts[r] - k < cols[1]) next
    # Move the k germline bases before the run to its tail; query columns stay put.
    window <- (starts[r] - k):ends[r]
    g[window] <- c(rep(DASH, len), g[(starts[r] - k):(starts[r] - 1)])
  }
  g
}

# Cuts a flank to whole codons at its far end, counting germline columns only.
trim_to_codons <- function(cols, g, from_start) {
  keep <- cols[g[cols] != DASH]
  d <- length(keep) %% 3
  if (!d) return(cols)
  if (from_start) cols[cols > keep[d]] else cols[cols < keep[length(keep) - d + 1]]
}

reframe <- function(seq, germ, junction) {
  s <- charToRaw(seq)
  g <- charToRaw(germ)
  if (length(s) != length(g) || !nzchar(junction)) return(NULL)
  # IMGT gaps out: inputs may be gapped to different references or not at all. A query base
  # under a germline gap is an insertion, dropped with the others below.
  g[g == DOT & s != DOT & s != DASH] <- DASH
  keep <- g != DOT
  gap_columns <- sum(!keep)
  s <- s[keep]
  g <- g[keep]
  bases <- which(s != DOT & s != DASH)
  if (!length(bases)) return(NULL)
  ungapped <- toupper(rawToChar(s[bases]))
  at <- regexpr(toupper(junction), ungapped, fixed = TRUE)
  # A mutated conserved codon can hide the junction; find the CDR3 and take a codon each side.
  if (at < 0 && nchar(junction) > 6) {
    core <- regexpr(toupper(substr(junction, 4, nchar(junction) - 3)), ungapped, fixed = TRUE)
    if (core >= 4 && core + nchar(junction) - 4 <= length(bases)) at <- core - 3
  }
  if (at < 0) return(NULL)
  from <- bases[at]
  to <- bases[at + nchar(junction) - 1]
  mid <- from:to
  mid <- mid[s[mid] != DASH]
  # A junction that is not whole codons would put the flanks in different frames.
  if (length(mid) %% 3 != 0) return(NULL)
  left <- if (from > bases[1]) bases[1]:(from - 1) else integer(0)
  right <- if (to < bases[length(bases)]) (to + 1):bases[length(bases)] else integer(0)
  left <- trim_to_codons(left, g, TRUE)
  right <- trim_to_codons(right, g, FALSE)
  # The slider moves insertions within a flank, never across its edge, so the frame holds.
  g <- codon_align_insertions(g, left)
  g <- codon_align_insertions(g, right)
  columns <- length(left) + length(right)
  left <- left[g[left] != DASH]
  right <- right[g[right] != DASH]
  n_mid <- rep(NB, length(mid))
  list(seq = rawToChar(c(s[left], s[mid], s[right])),
       germ = rawToChar(c(g[left], n_mid, g[right])),
       masked = rawToChar(c(s[left], n_mid, s[right])),
       dropped = columns - length(left) - length(right),
       left = length(left), right = length(right), gaps = gap_columns)
}

# In chunks, written straight into columns: one result list per row ran out of memory at 1M rows.
REFRAME_CHUNK <- as.integer(Sys.getenv("REFRAME_CHUNK", "50000"))
n <- nrow(joined)
seq_out <- germ_out <- masked_out <- rep(NA_character_, n)
left_out <- right_out <- dropped <- gaps <- rep(NA_integer_, n)
for (from in seq.int(1L, n, by = REFRAME_CHUNK)) {
  idx <- from:min(n, from + REFRAME_CHUNK - 1L)
  # Rows are independent, so a chunk is split over the stage's workers.
  # Contiguous blocks, so the results come back in row order.
  size <- ceiling(length(idx) / min(fork_workers(threads), length(idx)))
  parts <- split(idx, (seq_along(idx) - 1L) %/% size)
  framed <- unlist(parallel::mclapply(parts, function(rows) {
    mapply(reframe, joined$sequence_alignment[rows], joined$germline_alignment[rows],
           joined$junction[rows], SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }, mc.cores = length(parts)), recursive = FALSE, use.names = FALSE)
  ok <- idx[!vapply(framed, is.null, logical(1))]
  framed <- framed[!vapply(framed, is.null, logical(1))]
  pick <- function(field, type) vapply(framed, function(f) f[[field]], type)
  seq_out[ok] <- pick("seq", character(1))
  germ_out[ok] <- pick("germ", character(1))
  masked_out[ok] <- pick("masked", character(1))
  left_out[ok] <- pick("left", integer(1))
  right_out[ok] <- pick("right", integer(1))
  dropped[ok] <- pick("dropped", integer(1))
  gaps[ok] <- pick("gaps", integer(1))
}
placed <- !is.na(seq_out)
if (any(!placed)) {
  cat(sprintf("%d rows lost their alignment: the junction was not found in the aligned sequence, or was not whole codons\n",
              sum(!placed)))
}
joined <- joined[placed, , drop = FALSE]
if (!nrow(joined)) finish_empty("no aligned sequence contains its junction")
joined$sequence_alignment <- seq_out[placed]
joined$germline_alignment <- germ_out[placed]
# How far each row reaches on either side of its junction, for the padding below.
joined$frame_left <- left_out[placed]
joined$frame_right <- right_out[placed]
dropped <- dropped[placed]
gaps <- gaps[placed]
if (any(gaps > 0)) {
  cat(sprintf("IMGT gaps removed from %d of %d alignments (%d gap columns), so alignments gapped to different references, or not gapped, line up\n",
              sum(gaps > 0), length(gaps), sum(gaps)))
}
cat(sprintf("alignments rebuilt in the germline frame: %d insertion columns dropped from %d of %d rows\n",
            sum(dropped), sum(dropped > 0), length(dropped)))
# Coverage either side of the junction (full VDJRegion: 309 left); short rows are kept.
cat(sprintf("covered germline columns either side of the junction: left %d to %d, right %d to %d\n",
            min(joined$frame_left), max(joined$frame_left),
            min(joined$frame_right), max(joined$frame_right)))

# Junction masked in the query too, so HILARy's adaptive mode counts V and J mutations only.
joined$masked_sequence_alignment <- masked_out[placed]

write.table(joined[, ALIGNED_COLUMNS], aligned_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("aligned table written: %d rows, %d clonotypes\n", nrow(joined), length(unique(joined$sequence_id))))
progress("Alignments joined")
close_log()
quit(status = 0)
}

# Trees stage, from the align stage's table.
suppressMessages({library(dowser); library(ape); library(dplyr)})

if (is.null(aligned_in_path)) stop("missing --aligned")
joined <- read_tsv(aligned_in_path)
if (!nrow(joined)) finish_empty("no clonotype carries a heavy chain alignment")
joined$frame_left <- as.integer(joined$frame_left)
joined$frame_right <- as.integer(joined$frame_right)
if (use_light && !any(joined$locus != HEAVY)) {
  cat("no light chain alignments: light chain resolution off\n")
  use_light <- FALSE
}
cat(sprintf("read %d aligned rows for %d clonotypes\n", nrow(joined), length(unique(joined$sequence_id))))

# One tip per clonotype, so the clonotype key is the cell.
joined$cell_id <- joined$sequence_id
# resolveLightChains reads junction length from this column; without it the call fails.
joined$junction_length <- nchar(joined$junction)
joined$clone_id <- clones$clone_id[match(joined$sequence_id, clones$sequence_id)]
joined <- joined[!is.na(joined$clone_id), ]
if (!nrow(joined)) finish_empty("no aligned clonotype belongs to a clone")

# Lineage refinement.

if (use_light) {
  # minseq 1: a single-member clone still needs a lineage id.
  resolved <- resolveLightChains(joined, cell = "cell_id", locus = "locus",
                                 heavy = HEAVY, minseq = 1)
  resolved$lineage_id <- resolved$clone_subgroup_id
  h <- resolved[resolved$locus == HEAVY, ]
  cat(sprintf("resolveLightChains: %d clones into %d lineages, %d members without a light chain\n",
              length(unique(h$clone_id)), length(unique(h$lineage_id)),
              sum(h$vj_gene == "missing")))
} else {
  resolved <- joined
  resolved$lineage_id <- resolved$clone_id
}
# Dowser puts the clone id unquoted into temp file names and shell commands, and a
# donor name may hold "/" or a space. Build under flat surrogate ids, map back after.
lineage_ids <- unique(resolved$lineage_id)
surrogate <- setNames(paste0("L", seq_along(lineage_ids)), lineage_ids)
real_lineage <- setNames(lineage_ids, surrogate)
resolved$lineage_key <- surrogate[resolved$lineage_id]
# Clonotype ids get flat surrogates too: the tree files rewrite ":", ";", ",", "=" and spaces
# in tip names, so tips would no longer match. Heavy and light rows of one clonotype share one.
# tree_rows maps tips back; nothing else sees the surrogate.
tip_ids <- unique(resolved$sequence_id)
resolved$tip_id <- paste0("t", match(resolved$sequence_id, tip_ids))
TIP_REAL <- setNames(tip_ids, paste0("t", seq_along(tip_ids)))

# Subgroup is in lineage_id; flat, so formatClones neither appends it nor drops light rows.
resolved$clone_subgroup <- 1L

# Membership: every clustered clonotype gets a row; unmatched ones keep a lineage, no tree.

membership <- data.frame(
  sequence_id = clones$sequence_id,
  lineage_id = resolved$lineage_id[match(clones$sequence_id, resolved$sequence_id)],
  stringsAsFactors = FALSE
)
unmatched <- is.na(membership$lineage_id)
if (any(unmatched)) {
  # Subgroup 1 is the largest, where resolveLightChains puts unknown light chains.
  membership$lineage_id[unmatched] <- if (use_light) {
    paste0(clones$clone_id[unmatched], "_1")
  } else {
    clones$clone_id[unmatched]
  }
  cat(sprintf("%d of %d clonotypes (%.1f%%) had no AIRR row; placed in their clone's first lineage, no tree\n",
              sum(unmatched), length(unmatched), pct(sum(unmatched), length(unmatched))))
}
# Sequence groups: clonotypes identical on every chain are one group and one tip,
# compared on the rebuilt alignments. A group with an anchor is an anchor and is
# represented by it, else by the lowest id. Every member links to the tip.
group_of <- function(d) {
  chains <- split(paste(d$locus, d$sequence_alignment), d$sequence_id)
  signature <- vapply(chains, function(x) paste(sort(x), collapse = "|"), character(1))
  lineage <- vapply(split(d$lineage_id, d$sequence_id), function(x) x[1], character(1))
  ids <- names(signature)
  groups <- paste(lineage[ids], signature)
  # Anchors first, then lowest id; `order` is stable, so each group's first row represents it.
  o <- order(groups, !(ids %in% anchor_ids), ids)
  data.frame(sequence_id = ids[o], group_id = groups[o],
             is_representative = !duplicated(groups[o]),
             stringsAsFactors = FALSE)
}

progress("Grouping identical sequences")
groups <- group_of(resolved)
resolved$group_id <- groups$group_id[match(resolved$sequence_id, groups$sequence_id)]
resolved$is_representative <- groups$is_representative[match(resolved$sequence_id, groups$sequence_id)]
cat(sprintf("sequence groups: %d clonotypes are %d distinct sequences, %d of them seen more than once\n",
            length(unique(resolved$sequence_id)), length(unique(groups$group_id)),
            sum(table(groups$group_id) > 1)))

# The group lets collect size a lineage by distinct sequences (its tips).
membership$group_id <- resolved$group_id[match(membership$sequence_id, resolved$sequence_id)]
membership$group_id[is.na(membership$group_id)] <- membership$sequence_id[is.na(membership$group_id)]
# A linker p-column needs a value alongside its two key columns.
membership$link <- 1L
write.table(membership, lineages_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Heavy nt differences from the germline, V and J only, where both sides carry a base.
if (!is.null(germline_mutations_path)) {
  # Heavy rows first, then one per clonotype: a light row listed first must not win.
  h <- resolved[resolved$locus == HEAVY, , drop = FALSE]
  h <- h[!duplicated(h$sequence_id), , drop = FALSE]
  base <- c("A", "C", "G", "T")
  count_of <- function(q, g) {
    q <- strsplit(toupper(q), "")[[1]]
    g <- strsplit(toupper(g), "")[[1]]
    n <- seq_len(min(length(q), length(g)))
    sum(q[n] %in% base & g[n] %in% base & q[n] != g[n])
  }
  germline_mutations <- data.frame(
    sequence_id = h$sequence_id,
    lineage_id = h$lineage_id,
    germline_mutation_count = as.integer(mapply(count_of, h$masked_sequence_alignment,
                                                h$germline_alignment, USE.NAMES = FALSE)),
    stringsAsFactors = FALSE)
  write.table(germline_mutations, germline_mutations_path, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
  cat(sprintf("germline mutations: counted for %d clonotypes\n", nrow(germline_mutations)))
}

# aa-cdist: distance to the lineage's amino acid consensus (Ralph & Matsen 2020, PLoS
# Comput Biol 16(11):e1008391). Needs no tree, so uses every placed sequence. Heavy only.

# Below this the consensus is unreliable (two members score zero). Same floor as partis.
AA_CDIST_MIN_SEQUENCES <- 10

GAP_AA <- c(".", "-")
AMBIGUOUS_AA <- "X"

# Standard genetic code, codons in TCAG order on each of the three positions.
CODON_BASES <- c("T", "C", "A", "G")
AA_BY_CODON <- setNames(
  strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]],
  paste0(rep(CODON_BASES, each = 16), rep(rep(CODON_BASES, each = 4), 4),
         rep(CODON_BASES, 16)))

# Every codon's residue: pure gap is a gap; a non-base or a stop is X (never votes or scores).
CODON_TABLE <- {
  all_codons <- do.call(paste0, expand.grid(
    c(CODON_BASES, GAP_AA, AMBIGUOUS_AA), c(CODON_BASES, GAP_AA, AMBIGUOUS_AA),
    c(CODON_BASES, GAP_AA, AMBIGUOUS_AA), stringsAsFactors = FALSE)[, 3:1])
  aa <- vapply(all_codons, function(cd) {
    parts <- strsplit(cd, "")[[1]]
    if (all(parts %in% GAP_AA)) return(GAP_AA[1])
    if (any(!parts %in% CODON_BASES)) return(AMBIGUOUS_AA)
    a <- AA_BY_CODON[[cd]]
    if (a == "*") AMBIGUOUS_AA else a
  }, character(1), USE.NAMES = FALSE)
  setNames(aa, all_codons)
}
# The residues a branch can be said to have changed between; X and gaps never.
SETTLED_AA <- setdiff(unique(AA_BY_CODON), "*")

# Rows arrive rebuilt in their germline's frame, gap columns gone, so they translate in frame. Vectorised over rows.
translate_all <- function(nt) {
  nt <- toupper(nt)
  n <- nchar(nt) %/% 3
  width <- max(c(0L, n))
  out <- vector("list", length(nt))
  if (!width) return(lapply(seq_along(nt), function(i) character(0)))
  m <- matrix(GAP_AA[1], nrow = length(nt), ncol = width)
  for (k in seq_len(width)) {
    at <- n >= k
    if (!any(at)) next
    codons <- substring(nt[at], 3L * k - 2L, 3L * k)
    aa <- CODON_TABLE[codons]
    aa[is.na(aa)] <- AMBIGUOUS_AA
    m[at, k] <- aa
  }
  lapply(seq_along(nt), function(i) m[i, seq_len(n[i])])
}

# One vote per sequence; gaps and X do not vote, and a tie gives X, as in the paper.
aa_consensus <- function(m) {
  residues <- setdiff(unique(as.vector(m)), c(GAP_AA, AMBIGUOUS_AA))
  if (!length(residues)) return(rep(AMBIGUOUS_AA, ncol(m)))
  # One row per residue, one column per position, holding the votes cast.
  counts <- vapply(residues, function(r) colSums(m == r), numeric(ncol(m)))
  counts <- matrix(counts, nrow = ncol(m), dimnames = list(NULL, residues))
  best <- apply(counts, 1, max)
  winners <- rowSums(counts == best)
  out <- residues[max.col(counts, ties.method = "first")]
  out[best == 0 | winners > 1] <- AMBIGUOUS_AA
  out
}

# Score a position only where both residues are readable. Whole lineage at once.
aa_distances <- function(m, cons) {
  unreadable <- c(GAP_AA, AMBIGUOUS_AA)
  keep_col <- !cons %in% unreadable
  if (!any(keep_col)) return(rep(0L, nrow(m)))
  sub <- m[, keep_col, drop = FALSE]
  ref <- rep(cons[keep_col], each = nrow(sub))
  as.integer(rowSums(sub != ref & !(sub %in% unreadable)))
}

progress("Amino acid consensus distance")
heavy_seqs <- resolved[resolved$locus == HEAVY,
                       c("sequence_id", "lineage_id", "sequence_alignment", "frame_left", "group_id")]
heavy_seqs <- heavy_seqs[!duplicated(heavy_seqs$sequence_id), ]
# One vote per distinct sequence, as the paper does not weight by count.
voters <- heavy_seqs[!duplicated(heavy_seqs$group_id), ]

# Row indices per lineage in first-appearance order: one pass instead of a scan per lineage.
rows_by_lineage <- function(d) {
  split(seq_len(nrow(d)), factor(d$lineage_id, levels = unique(d$lineage_id)))
}

heavy_lineages <- rows_by_lineage(voters)
# The floor counts voters, not clonotypes.
support <- data.frame(lineage_id = as.character(names(heavy_lineages)),
                      consensus_sequence_count = unname(lengths(heavy_lineages)),
                      stringsAsFactors = FALSE)
members_by_group <- split(heavy_seqs$sequence_id, heavy_seqs$group_id)
# Plain vectors: slicing the data frame per lineage cost more than the work.
voter_seq <- voters$sequence_alignment
voter_frame <- as.integer(voters$frame_left)
voter_group <- voters$group_id

# Matched once: indexing a large named list by name rehashes on every call.
member_of <- match(voter_group, names(members_by_group))
scored_k <- which(lengths(heavy_lineages) >= AA_CDIST_MIN_SEQUENCES)
TRANSLATE_BLOCK <- 100000L
aa_of <- vector("list", nrow(voters))
scored_rows <- as.integer(unlist(heavy_lineages[scored_k], use.names = FALSE))
for (b in split(scored_rows, (seq_along(scored_rows) - 1L) %/% TRANSLATE_BLOCK)) {
  aa_of[b] <- translate_all(voter_seq[b])
}
rm(scored_rows)

cdist_parts <- vector("list", length(scored_k))
for (j in seq_along(scored_k)) {
  k <- scored_k[j]
  lid <- names(heavy_lineages)[k]
  rows <- heavy_lineages[[k]]
  aa <- aa_of[rows]
  # Align members on the junction: pad by the codons their 5' side is short.
  frame <- voter_frame[rows]
  lead <- (max(frame) - frame) %/% 3
  width <- max(lengths(aa) + lead)
  m <- matrix(GAP_AA[1], nrow = length(aa), ncol = width)
  for (i in seq_along(aa)) {
    if (length(aa[[i]])) m[i, lead[i] + seq_along(aa[[i]])] <- aa[[i]]
  }
  cons <- aa_consensus(m)
  scored <- aa_distances(m, cons)
  members <- members_by_group[member_of[rows]]
  cdist_parts[[j]] <- data.frame(
    sequence_id = unlist(members, use.names = FALSE),
    lineage_id = lid,
    aa_cdist = rep(scored, lengths(members)),
    stringsAsFactors = FALSE)
}

cdist <- if (length(cdist_parts)) as.data.frame(dplyr::bind_rows(cdist_parts)) else
  data.frame(sequence_id = character(0), lineage_id = character(0), aa_cdist = integer(0))
cat(sprintf("aa-cdist: %d of %d lineages reached %d distinct sequences, covering %d clonotypes\n",
            length(cdist_parts), nrow(support), AA_CDIST_MIN_SEQUENCES, nrow(cdist)))

if (!is.null(support_path)) {
  write.table(support, support_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
if (!is.null(cdist_path)) {
  write.table(cdist, cdist_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# Per-lineage downsampling, off by default (FastTree+RAxML is near-linear in tips).
# Random, not by abundance, which would oversample the expanded clade; counts distinct
# sequences, keeps anchors, seeded by lineage id. Dropped clonotypes lose only tree columns.

lineage_seed <- function(lid) sum(utf8ToInt(lid) * seq_along(utf8ToInt(lid))) %% .Machine$integer.max
progress("Selecting tips and building lineage germlines")
capped <- 0L
dropped_tips <- 0L
kept_rows <- rep(TRUE, nrow(resolved))
lineages <- rows_by_lineage(resolved)
for (k in seq_along(lineages)) {
  rows <- lineages[[k]]
  heavy <- rows[resolved$locus[rows] == HEAVY]
  heavy <- heavy[!duplicated(resolved$sequence_id[heavy])]
  # A group is a tip, so the cap counts groups.
  by_group <- split(resolved$sequence_id[heavy], resolved$group_id[heavy])
  if (length(by_group) <= max_tips) next
  with_anchor <- vapply(split(resolved$sequence_id[heavy] %in% anchor_ids, resolved$group_id[heavy]), any, logical(1))
  others <- which(!with_anchor)
  set.seed(lineage_seed(names(lineages)[k]))
  chosen <- c(which(with_anchor), others[sample.int(length(others), max(0L, max_tips - sum(with_anchor)))])
  keep <- unlist(by_group[chosen], use.names = FALSE)
  kept_rows[rows[!(resolved$sequence_id[rows] %in% keep)]] <- FALSE
  capped <- capped + 1L
  dropped_tips <- dropped_tips + length(by_group) - length(chosen)
}
if (capped > 0) {
  cat(sprintf("downsampling: %d lineages above %d distinct sequences subsampled, %d of %d distinct sequences (%.1f%%) left out of their trees\n",
              capped, max_tips, dropped_tips, length(unique(groups$group_id)),
              pct(dropped_tips, length(unique(groups$group_id)))))
  resolved <- resolved[kept_rows, , drop = FALSE]
}

# Germlines, per lineage and locus: dowser needs one germline per clone per locus.

# Consensus germline: disagreeing columns become N; N (padding) abstains.
consensus_germline <- function(g) {
  u <- unique(g)
  if (length(u) == 1) return(u)
  m <- do.call(rbind, strsplit(g, ""))
  paste(apply(m, 2, function(col) {
    votes <- unique(col[col != "N"])
    if (length(votes) == 1) votes else "N"
  }), collapse = "")
}

STOPS <- c("TAA", "TAG", "TGA")
codon_starts <- function(g, frame) seq(frame + 1, length.out = max(0, (nchar(g) - frame) %/% 3), by = 3)
stop_positions <- function(g, frame = 0) {
  i <- codon_starts(g, frame)
  i[toupper(substring(g, i, i + 2L)) %in% STOPS]
}

# igphyml refuses in-frame stops; mask any left in the germline (pseudogene, bad frame).
mask_stops <- function(g, label) {
  hits <- stop_positions(g, 0)
  if (!length(hits)) return(g)
  s <- strsplit(g, "")[[1]]
  for (k in hits) s[k:(k + 2)] <- "N"
  cat(sprintf("%s: masked %d germline stop codon(s)\n", label, length(hits)))
  paste(s, collapse = "")
}

frames <- list()
lineages <- rows_by_lineage(resolved)
skip_flag <- logical(length(lineages))
r_heavy <- resolved$locus == HEAVY
r_group <- resolved$group_id
for (k in seq_along(lineages)) {
  lid <- names(lineages)[k]
  rows <- lineages[[k]]
  # Needs two distinct sequences; one sequence is one tip and no tree.
  if (length(unique(r_group[rows][r_heavy[rows]])) < 2) { skip_flag[k] <- TRUE; next }
  lin <- resolved[rows, ]
  parts <- list()
  ok <- TRUE
  for (loc in unique(lin$locus)) {
    d <- lin[lin$locus == loc, ]
    # Pad each side with N to the longest member so germline positions line up.
    pad <- function(x, before, after) {
      paste0(strrep("N", before), x, strrep("N", after))
    }
    left_pad <- max(d$frame_left) - d$frame_left
    right_pad <- max(d$frame_right) - d$frame_right
    if (any(left_pad > 0 | right_pad > 0)) {
      d$sequence_alignment <- pad(d$sequence_alignment, left_pad, right_pad)
      d$germline_alignment <- pad(d$germline_alignment, left_pad, right_pad)
      cat(sprintf("lineage %s: %s members padded to a common frame (%d rows)\n",
                  lid, loc, sum(left_pad > 0 | right_pad > 0)))
    }
    if (length(unique(nchar(d$germline_alignment))) != 1 ||
        length(unique(nchar(d$sequence_alignment))) != 1) {
      cat(sprintf("lineage %s skipped: %s alignment lengths differ\n", lid, loc))
      ok <- FALSE
      break
    }
    g <- consensus_germline(d$germline_alignment)
    d$germline_alignment_d_mask <- mask_stops(g, paste(lid, loc))
    parts[[loc]] <- d
  }
  if (!ok) { skip_flag[k] <- TRUE; next }
  frames[[length(frames) + 1]] <- if (length(parts) == 1) parts[[1]] else do.call(rbind, parts)
}
skipped <- names(lineages)[skip_flag]
rm(r_heavy, r_group, skip_flag)
if (!length(frames)) finish_empty("no lineage had two or more members with usable alignments",
                                  membership)
db <- as.data.frame(dplyr::bind_rows(frames))

# Trees.

chain <- if (use_light && any(db$locus != HEAVY)) "HL" else "H"
# `collapse = FALSE`: groups were deduplicated above, so dowser must not pick which id survives.
progress("Preparing lineages for tree building")
tips_db <- db[db$is_representative, , drop = FALSE]

# formatClones runs per lineage in the pool, and once here for IgPhyML lineages.
format_lineages <- function(d, nproc) {
  formatClones(d, clone = "lineage_key", seq = "sequence_alignment",
               germ = "germline_alignment_d_mask", id = "tip_id",
               cell = "tip_id", locus = "locus", heavy = HEAVY,
               chain = chain, split_light = FALSE, minseq = 2,
               collapse = FALSE, nproc = nproc)
}

# Both chains: per-chain omega and rate (IgPhyML), scaled branch lengths (RAxML).
partition_for <- function(b) if (chain == "HL") { if (b == "igphyml") "hl" else "scaled" } else NULL

# Builder per lineage: IgPhyML (slow) for none, anchor lineages, or all; RAxML otherwise.
tip_rows <- split(seq_len(nrow(tips_db)), factor(tips_db$lineage_key, levels = unique(tips_db$lineage_key)))
lineage_keys <- names(tip_rows)
anchored_keys <- unique(resolved$lineage_key[resolved$sequence_id %in% anchor_ids])
builder_of <- if (igphyml_scope == "all") {
  rep("igphyml", length(lineage_keys))
} else if (igphyml_scope == "anchored") {
  ifelse(lineage_keys %in% anchored_keys, "igphyml", "raxml")
} else {
  rep("raxml", length(lineage_keys))
}
cat(sprintf("tree builders: %d lineages with IgPhyML, %d with FastTree+RAxML (scope %s, %d lineages hold an anchor)\n",
            sum(builder_of == "igphyml"), sum(builder_of == "raxml"), igphyml_scope,
            sum(lineage_keys %in% anchored_keys)))

# Ancestral sequences: IgPhyML and parsimony write IUPAC codes where uncertain, dowser's
# RAxML the single best base. So RAxML gets IgPhyML's credible-set rule, applied to
# raxml-ng's `.ancestralProbs`. Parsimony keeps its flat threshold (not probabilities).

# Result names per builder; only the branch length unit differs.
BUILDER_LABEL <- c(raxml = "fasttree-raxml", igphyml = "igphyml", pratchet = "pratchet")

ASR_CREDIBLE_MASS <- 0.95

# IUPAC code per subset of {A,C,G,T}, indexed by a bitmask A=1 C=2 G=4 T=8.
IUPAC_BY_MASK <- c("N", "A", "C", "M", "G", "R", "S", "V",
                   "T", "W", "Y", "H", "K", "D", "B", "N")

# Settled bases; anything else is an ambiguity code and counts as unsettled.
UNAMBIGUOUS <- c("A", "C", "G", "T")

# Per-site marginals to one IUPAC string; identical rows are handled once.
iupac_from_probs <- function(p) {
  key <- paste(p[, 1], p[, 2], p[, 3], p[, 4])
  first <- !duplicated(key)
  u <- p[first, , drop = FALSE]
  sorted <- t(apply(u, 1, sort, decreasing = TRUE))
  cumulative <- t(apply(sorted, 1, cumsum))
  reached <- cumulative >= ASR_CREDIBLE_MASS
  # The full set always qualifies, so rounding short of the mass gives full ambiguity.
  reached[, ncol(reached)] <- TRUE
  crossed <- max.col(reached, ties.method = "first")
  cutoff <- sorted[cbind(seq_len(nrow(sorted)), crossed)]
  keep <- u >= cutoff
  mask <- keep[, 1] + 2 * keep[, 2] + 4 * keep[, 3] + 8 * keep[, 4]
  paste(IUPAC_BY_MASK[mask + 1][match(key, key[first])], collapse = "")
}

# Swap RAxML node sequences for credible-set ones; a missing or bad file keeps dowser's.
rethreshold_raxml <- function(built, dir, run_id) {
  rewritten <- 0
  for (i in seq_len(nrow(built))) {
    clone <- as.character(built$clone_id[i])
    prefix <- paste0(run_id, "_", clone, "_")
    path <- file.path(dir, paste0(prefix, "1_asr.raxml.ancestralProbs"))
    done <- tryCatch({
      if (!file.exists(path)) stop("no probability file")
      # Partitioned: a Part column, sites numbered per partition; part then site is alignment order.
      partitioned <- "Part" %in% trimws(strsplit(readLines(path, n = 1), "\t")[[1]])
      probs <- read.delim(path, sep = "\t", colClasses = c(
        "character", if (partitioned) "integer", "integer", "NULL",
        "numeric", "numeric", "numeric", "numeric"))
      if (!partitioned) probs$Part <- 1L
      tree <- built$trees[[i]]
      internal <- seq.int(length(tree$tip.label) + 1, length(tree$nodes))
      node_label <- vapply(internal, function(n) {
        named <- names(tree$nodes[[n]]$sequence)
        if (is.null(named)) NA_character_ else sub("^sequence\\.", "", named[1])
      }, character(1))
      for (node in split(seq_len(nrow(probs)), probs$Node)) {
        hit <- internal[which(node_label == probs$Node[node[1]])]
        if (!length(hit)) next
        rows <- node[order(probs$Part[node], probs$Site[node])]
        # RAxML writes "-" where the alignment had a gap; Dowser calls that N and so do we.
        replacement <- gsub("-", "N", iupac_from_probs(
          as.matrix(probs[rows, c("p_A", "p_C", "p_G", "p_T"), drop = FALSE])))
        if (nchar(replacement) != nchar(as.character(tree$nodes[[hit[1]]]$sequence)[1])) {
          stop("probability sites do not cover the node sequence")
        }
        for (n in hit) {
          names(replacement) <- names(tree$nodes[[n]]$sequence)
          tree$nodes[[n]]$sequence <- replacement
        }
      }
      built$trees[[i]] <- tree
      rm(probs)
      TRUE
    }, error = function(e) {
      cat(sprintf("%s: keeping RAxML's own node sequences (%s)\n", clone, conditionMessage(e)))
      FALSE
    })
    if (done) rewritten <- rewritten + 1
    unlink(list.files(dir, pattern = paste0("^", prefix), full.names = TRUE))
  }
  # Report only lineages that could not be thresholded; each prints its reason above.
  if (rewritten < nrow(built)) {
    cat(sprintf("ancestral sequences kept unthresholded for %d of %d lineages\n",
                nrow(built) - rewritten, nrow(built)))
  }
  built
}

# Tree pool: fresh R sessions, not forks (a fork's GC touches the whole parent heap),
# fed lineages as they free up. Workers return node rows, not trees, to save memory.

# IgPhyML and RAxML pick a seed from the clock unless given one; Dowser defaults RAxML to 28.
TREE_SEED <- 0L

attempt_build <- function(build, p, on, nproc) {
  if (!nrow(on)) return(simpleError("no eligible lineages"))
  tryCatch(
    if (build == "igphyml") {
      # Also buildIgphyml's default; named to show the shared threshold.
      getTrees(on, build = "igphyml", exec = igphyml, nproc = nproc, quiet = 1,
               partition = p, asrc = ASR_CREDIBLE_MASS, rseed = TREE_SEED)
    } else if (build == "raxml") {
      # Dowser deletes RAxML's working files, which hold the marginals, so use our own dir.
      run_id <- "lt"
      dir <- tempfile("raxml-")
      on.exit(unlink(dir, recursive = TRUE), add = TRUE)
      rethreshold_raxml(
        getTrees(on, build = "raxml", exec = raxml, nproc = nproc, quiet = 1,
                 partition = p, dir = dir, id = run_id, rm_temp = FALSE, rseed = TREE_SEED),
        dir, run_id)
    } else {
      getTrees(on, build = "pratchet", nproc = nproc, quiet = 1)
    },
    error = function(e) e)
}
usable_build <- function(r) !inherits(r, "error") && nrow(r) > 0
why_build <- function(r) if (inherits(r, "error")) conditionMessage(r) else "no trees returned"

# One lineage: the builder, unpartitioned if that fails, parsimony if that fails too.
build_lineage <- function(one, build, p) {
  r <- attempt_build(build, p, one, 1L)
  if (usable_build(r)) return(list(tree = r, builder = build))
  first <- why_build(r)
  if (!is.null(p)) {
    r <- attempt_build(build, NULL, one, 1L)
    if (usable_build(r)) return(list(tree = r, builder = build, note = sprintf("rebuilt unpartitioned (%s)", first)))
  }
  if (build != "pratchet") {
    r <- attempt_build("pratchet", NULL, one, 1L)
    if (usable_build(r)) {
      return(list(tree = r, builder = "pratchet",
                  note = sprintf("built with maximum parsimony rather than %s (%s)", build, first)))
    }
  }
  list(tree = NULL, note = sprintf("produced no tree (%s)", why_build(r)))
}

# Merge chains of internal nodes with identical sequences: the builder forced the
# split, no ancestor was inferred. Tips are kept. A failed collapse keeps the tree.
collapse_tree <- function(p) {
  collapsed <- tryCatch(collapseNodes(p), error = function(e) e)
  if (!inherits(collapsed, "error")) return(collapsed)
  cat(sprintf("%s: node collapse failed (%s); keeping the uncollapsed tree\n",
              p$name, conditionMessage(collapsed)))
  p
}

# Node rows and anchor distances for one tree. `lid` is the real lineage id, `gapped_*`
# the rebuilt rows with their deletions, `heavy_width` where heavy ends (dowser joins heavy then light).
tree_rows <- function(p, lid, gapped_tips, gapped_germ, heavy_width = NA_integer_) {
  anchors <- NULL
  # Tips carry surrogate ids (see TIP_REAL); "Germline" and anything unmapped stay as they are.
  tips <- p$tip.label
  real <- unname(TIP_REAL[tips])
  labels <- c(ifelse(is.na(real), tips, real), rep(NA_character_, p$Nnode))
  parent <- rep(NA_integer_, length(labels))
  dist <- rep(NA_real_, length(labels))
  # A collapsed tree's edge matrix can come back double; the anchor walk wants ids.
  parent[p$edge[, 2]] <- as.integer(p$edge[, 1])
  dist[p$edge[, 2]] <- p$edge.length
  # Reconstructed sequence at each node; a tip's is its observed one.
  sequences <- vapply(seq_along(labels), function(n) {
    s <- if (n <= length(p$nodes)) p$nodes[[n]]$sequence else NULL
    if (is.null(s) || !length(s)) "" else unname(as.character(s)[1])
  }, character(1))
  # Tips and the germline take back their gapped originals where the frame matches;
  # inferred nodes get a gap only where every tip has one, N stays for the rest.
  for (n in which(!is.na(labels))) {
    orig <- if (labels[n] == "Germline") gapped_germ else gapped_tips[labels[n]]
    if (!is.na(orig) && nchar(orig) == nchar(sequences[n])) sequences[n] <- unname(orig)
  }
  tips <- which(!is.na(labels) & labels != "Germline" & nzchar(sequences))
  if (length(tips) && length(unique(nchar(sequences[tips]))) == 1) {
    tip_bases <- do.call(rbind, strsplit(sequences[tips], ""))
    all_gap <- apply(tip_bases, 2, function(col) all(col %in% c(".", "-")))
    for (n in which(is.na(labels) & nchar(sequences) == ncol(tip_bases))) {
      s <- strsplit(sequences[n], "")[[1]]
      s[all_gap] <- tip_bases[1, all_gap]
      sequences[n] <- paste(s, collapse = "")
    }
  }
  # Depth and root distance per node; edges parent first, so one pass does it.
  ordered <- reorder(p, "cladewise")
  edge <- ordered$edge
  edge_length <- if (is.null(ordered$edge.length)) rep(NA_real_, nrow(edge)) else ordered$edge.length
  depth <- integer(length(labels))
  root_dist <- numeric(length(labels))
  for (e in seq_len(nrow(edge))) {
    depth[edge[e, 2]] <- depth[edge[e, 1]] + 1L
    root_dist[edge[e, 2]] <- root_dist[edge[e, 1]] + edge_length[e]
  }
  # Observed tips below each node, children before parents; the germline is not one.
  below <- as.integer(!is.na(labels) & labels != "Germline")
  post <- reorder(p, "postorder")$edge
  for (e in seq_len(nrow(post))) below[post[e, 1]] <- below[post[e, 1]] + below[post[e, 2]]
  # Own branch over root distance: a unitless ratio, comparable across builders.
  terminal_fraction <- ifelse(!is.na(dist) & !is.na(root_dist) & root_dist > 0, dist / root_dist, NA_real_)
  parent_descendants <- ifelse(is.na(parent), NA_integer_, below[parent])

  # A mutation needs a settled base at both ends; otherwise it counts as unresolved.
  # Counted per chain (1-based in the chain's rebuilt frame, not IMGT) and per alphabet;
  # a chain that is not whole codons gets no amino acid figures.
  width <- max(c(0L, nchar(sequences)))
  heavy_end <- if (!is.na(heavy_width) && heavy_width > 0L && heavy_width < width) heavy_width else width
  chain_spans <- list(heavy = c(1L, heavy_end))
  if (width > heavy_end) chain_spans$light <- c(heavy_end + 1L, width)

  n_nodes <- length(labels)
  # Settled changes on a branch: both ends in `alphabet`. Other differences are unresolved.
  branch_step <- function(from, to, alphabet) {
    if (!length(from) || length(from) != length(to)) return(list(text = "", count = 0L, unresolved = 0L))
    differs <- which(from != to)
    settled <- from[differs] %in% alphabet & to[differs] %in% alphabet
    list(text = paste(paste0(from[differs][settled], differs[settled], to[differs][settled]),
                      collapse = ","),
         count = sum(settled), unresolved = sum(!settled))
  }
  steps <- list()
  # Whether a node carries the chain: a missing light chain is all N.
  carries <- list()
  for (chain in names(chain_spans)) {
    span <- chain_spans[[chain]]
    part <- toupper(substr(sequences, span[1], span[2]))
    part[nchar(sequences) != width] <- ""
    nt <- lapply(part, function(x) if (nzchar(x)) strsplit(x, "")[[1]] else character(0))
    in_frame <- (span[2] - span[1] + 1L) %% 3L == 0L
    if (!in_frame) {
      cat(sprintf("%s: %s chain is %d columns, not whole codons; no amino acid figures for it\n",
                  lid, chain, span[2] - span[1] + 1L))
    }
    aa <- if (in_frame) translate_all(part) else NULL
    s <- list(nt_text = rep("", n_nodes), nt_count = rep(0L, n_nodes),
              unresolved = rep(0L, n_nodes),
              aa_text = rep(if (in_frame) "" else NA_character_, n_nodes),
              aa_count = rep(if (in_frame) 0L else NA_integer_, n_nodes))
    for (n in seq_len(n_nodes)) {
      up <- parent[n]
      if (is.na(up)) next
      b <- branch_step(nt[[up]], nt[[n]], UNAMBIGUOUS)
      s$nt_text[n] <- b$text
      s$nt_count[n] <- b$count
      s$unresolved[n] <- b$unresolved
      if (in_frame) {
        b <- branch_step(aa[[up]], aa[[n]], SETTLED_AA)
        s$aa_text[n] <- b$text
        s$aa_count[n] <- b$count
      }
    }
    carried <- vapply(nt, function(x) any(x %in% UNAMBIGUOUS), logical(1))
    # A light half no node settles is padding, not a chain: blank, not zero.
    if (chain == "light" && !any(carried)) next
    steps[[chain]] <- s
    carries[[chain]] <- carried
  }

  # Anchored search: settled mutations along the tree from each non-anchor tip to the
  # nearest anchor (by heavy aa, then nt, then label). Light figures need light at both ends.
  anchors_here <- which(labels %in% anchor_ids)
  if (length(anchors_here)) {
    # Settled mutations from the root to each node, so any path cost is a difference.
    root_cost <- function(per_branch) {
      cost <- rep(0L, n_nodes)
      for (e in seq_len(nrow(edge))) cost[edge[e, 2]] <- cost[edge[e, 1]] + per_branch[edge[e, 2]]
      cost
    }
    costs <- lapply(steps, function(s) list(nt = root_cost(s$nt_count), aa = root_cost(s$aa_count)))
    on_line <- function(n) {
      out <- logical(n_nodes)
      while (!is.na(n)) { out[n] <- TRUE; n <- parent[n] }
      out
    }
    anchors_here <- anchors_here[order(labels[anchors_here])]
    anchor_lines <- lapply(anchors_here, on_line)
    candidates <- which(!is.na(labels) & labels != "Germline" & !(labels %in% anchor_ids))
    if (length(candidates)) {
      path <- function(cost, cand, k, lca) cost[cand] + cost[anchors_here[k]] - 2L * cost[lca]
      rows <- lapply(candidates, function(cand) {
        lcas <- vapply(seq_along(anchors_here), function(k) {
          lca <- cand
          while (!anchor_lines[[k]][lca]) lca <- parent[lca]
          lca
        }, integer(1))
        heavy_nt <- vapply(seq_along(anchors_here), function(k) path(costs$heavy$nt, cand, k, lcas[k]), integer(1))
        heavy_aa <- vapply(seq_along(anchors_here), function(k) path(costs$heavy$aa, cand, k, lcas[k]), integer(1))
        # `order` is stable, so ties fall to the label order set above.
        k <- order(heavy_aa, heavy_nt)[1]
        light_ok <- !is.null(costs$light) && carries$light[cand] && carries$light[anchors_here[k]]
        data.frame(
          sequence_id = labels[cand],
          anchor_id = labels[anchors_here[k]],
          anchor_aa_heavy = heavy_aa[k],
          anchor_nt_heavy = heavy_nt[k],
          anchor_aa_light = if (light_ok) path(costs$light$aa, cand, k, lcas[k]) else NA_integer_,
          anchor_nt_light = if (light_ok) path(costs$light$nt, cand, k, lcas[k]) else NA_integer_,
          stringsAsFactors = FALSE)
      })
      anchors <- do.call(rbind, rows)
    }
  }

  nodes <- data.frame(
    lineage_id = lid,
    node_id = seq_along(labels),
    parent_id = parent,
    distance = dist,
    is_observed = ifelse(is.na(labels) | labels == "Germline", "false", "true"),
    label = ifelse(is.na(labels), "", labels),
    heavy_sequence = substr(sequences, 1L, heavy_end),
    light_sequence = if (is.null(steps$light)) NA_character_ else substr(sequences, heavy_end + 1L, width),
    node_depth = depth,
    terminal_branch_fraction = terminal_fraction,
    parent_descendant_count = parent_descendants,
    stringsAsFactors = FALSE
  )
  for (chain in c("heavy", "light")) {
    s <- steps[[chain]]
    blank_text <- rep(NA_character_, n_nodes)
    blank_count <- rep(NA_integer_, n_nodes)
    nodes[[paste0(chain, "_mutations_from_parent")]] <- if (is.null(s)) blank_text else s$nt_text
    nodes[[paste0(chain, "_mutation_count_from_parent")]] <- if (is.null(s)) blank_count else s$nt_count
    nodes[[paste0(chain, "_unresolved_from_parent")]] <- if (is.null(s)) blank_count else s$unresolved
    nodes[[paste0(chain, "_aa_mutations_from_parent")]] <- if (is.null(s)) blank_text else s$aa_text
    nodes[[paste0(chain, "_aa_mutation_count_from_parent")]] <- if (is.null(s)) blank_count else s$aa_count
  }
  list(nodes = nodes[, NODE_COLUMNS], anchors = anchors)
}

# A build step to its result: collapsed tree rows and builder.
finish_lineage <- function(step, u) {
  out <- c(list(slot = u$slot, label = if (is.null(step$label)) u$label else step$label,
                note = step$note), step$meta)
  if (is.null(step$tree)) return(out)
  p <- step$tree$trees[[1]]
  before <- p$Nnode
  p <- collapse_tree(p)
  # dowser's record says where the heavy chain ends; our gapped germline is the fallback.
  heavy_width <- tryCatch(nchar(step$tree$data[[1]]@germline), error = function(e) NA_integer_)
  if (is.na(heavy_width) && !is.null(u$gapped_germ) && !is.na(u$gapped_germ)) heavy_width <- nchar(u$gapped_germ)
  c(out, tree_rows(p, u$lid, u$gapped_tips, u$gapped_germ, heavy_width),
    list(clone_id = as.character(step$tree$clone_id[1]),
         builder = BUILDER_LABEL[[step$builder]],
         nodes_before = before, nodes_after = p$Nnode))
}

# Formats one lineage and picks its build. Under three tips goes to parsimony: FastTree
# and raxml-ng need four sequences with the germline, and three taxa have one topology.
route_lineage <- function(u) {
  fmt <- tryCatch(format_lineages(u$db, 1L), error = function(e) e)
  if (inherits(fmt, "error") || !nrow(fmt)) {
    # A lineage formatClones empties is dropped, as the single call dropped it.
    if (inherits(fmt, "error") && !grepl("No clones remain", conditionMessage(fmt))) {
      cat(sprintf("%s: formatClones failed (%s)\n", u$lid, conditionMessage(fmt)))
    }
    return(list(tree = NULL, meta = list(formatted = FALSE)))
  }
  one <- fmt[1, ]
  mixed <- chain == "HL" && any(one$data[[1]]@locus != HEAVY)
  meta <- list(formatted = TRUE, kept_seqs = as.integer(one$seqs), mixed = mixed,
               below_min = min_tips > 2 && one$seqs < min_tips, small = FALSE)
  label <- paste("raxml", if (mixed) "with light chains" else "heavy chain only")
  if (meta$below_min) return(list(tree = NULL, meta = meta))
  build <- "raxml"
  part <- if (mixed) partition_for("raxml") else NULL
  if (one$seqs < 3) {
    build <- "pratchet"
    part <- NULL
    meta$small <- TRUE
  }
  step <- build_lineage(one, build, part)
  step$label <- label
  step$meta <- meta
  step
}

# One lineage's worker job; output, warnings and messages are captured and returned.
process_unit <- function(u) {
  # Parsimony and its random resolution draw from R's RNG; a worker's own seed is the clock's.
  set.seed(lineage_seed(u$lid))
  out <- NULL
  # formatClones warns per call about stop-codon removals; sum the counts instead.
  stops <- 0L
  log <- utils::capture.output(out <- withCallingHandlers(
    finish_lineage(if (is.null(u$db)) build_lineage(u$one, u$build, u$part) else route_lineage(u), u),
    warning = function(w) {
      m <- conditionMessage(w)
      hit <- regmatches(m, regexec("^([0-9]+) sequence\\(s\\) with an inframe stop codon were removed", m))[[1]]
      if (length(hit)) stops <<- stops + as.integer(hit[2]) else cat(sprintf("warning: %s\n", m))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cat(conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  out$log <- log
  out$stops_removed <- stops
  out
}

# What a worker needs besides dowser: the functions above and what they read.
POOL_EXPORTS <- c("attempt_build", "usable_build", "why_build", "build_lineage",
                  "collapse_tree", "tree_rows", "finish_lineage", "process_unit",
                  "route_lineage", "format_lineages", "partition_for", "chain", "HEAVY",
                  "lineage_seed", "TREE_SEED", "TIP_REAL",
                  "min_tips",
                  "rethreshold_raxml", "iupac_from_probs", "igphyml", "raxml",
                  "ASR_CREDIBLE_MASS", "IUPAC_BY_MASK", "UNAMBIGUOUS", "BUILDER_LABEL",
                  "anchor_ids",
                  # tree_rows counts steps per chain and in amino acids.
                  "NODE_COLUMNS", "SETTLED_AA", "translate_all", "CODON_TABLE", "GAP_AA",
                  "AMBIGUOUS_AA")

# Workers settle near 1.9 GiB at any lineage size and reach 2.8 GiB on the largest;
# 1.5 GiB let the pool fill 83-97% of its limit.
WORKER_BUDGET <- 2.5 * 2^30
fit_pool <- function(want) {
  limit <- memory_limit()
  rss <- self_rss()
  if (is.na(limit) || is.na(rss)) {
    cat(sprintf("tree builders: %d workers (no memory limit readable, so not bounded)\n", want))
    return(want)
  }
  fits <- max(1L, as.integer(floor((limit * 0.8 - rss) / WORKER_BUDGET)))
  got <- min(want, fits)
  cat(sprintf("tree builders: %d of %d workers (%.1f GiB limit, %.1f GiB held here, %.1f GiB budgeted per worker)\n",
              got, want, limit / 2^30, rss / 2^30, WORKER_BUDGET / 2^30))
  got
}

# Largest lineage first, so a big one is not the run's tail; clusterApplyLB's loop, inlined.
run_pool <- function(cl, units) {
  n <- length(units)
  results <- vector("list", n)
  queue <- order(-vapply(units, function(u) as.numeric(u$size), numeric(1)))
  sent <- 0L
  send <- function(node) {
    sent <<- sent + 1L
    parallel:::sendCall(cl[[node]], process_unit, list(units[[queue[sent]]]), tag = queue[sent])
  }
  for (node in seq_len(min(length(cl), n))) send(node)
  for (got in seq_len(n)) {
    r <- parallel:::recvOneResult(cl)
    if (sent < n) send(r$node)
    value <- r$value
    if (inherits(value, "try-error")) {
      value <- list(slot = units[[r$tag]]$slot, label = units[[r$tag]]$label,
                    note = sprintf("produced no tree (%s)", trimws(as.character(value))))
      # A lineage routed in the worker has no label yet, so its note would be lost.
      if (is.null(value$label)) cat(sprintf("%s: %s\n", units[[r$tag]]$lid, value$note))
    }
    if (length(value$log)) cat(value$log, sep = "\n")
    value$log <- NULL
    results[[r$tag]] <- value
    trees_done <<- trees_done + units[[r$tag]]$size
    now <- Sys.time()
    if (as.numeric(difftime(now, last_tree_progress, units = "secs")) >= TREE_PROGRESS_EVERY || got == n) {
      last_tree_progress <<- now
      progress(sprintf("Trees: %.1f%%", 100 * trees_done / max(1, trees_total)))
    }
  }
  results
}

# Free what builders do not need; keep gapped heavy originals (dowser masks gaps to N).
heavy_tips <- tips_db[tips_db$locus == HEAVY & !duplicated(tips_db$sequence_id), ]
tips_by_lineage <- split(setNames(toupper(heavy_tips$sequence_alignment), heavy_tips$sequence_id),
                         heavy_tips$lineage_key)
gapped_germline <- with(heavy_tips[!duplicated(heavy_tips$lineage_key), ],
                        setNames(toupper(germline_alignment_d_mask), lineage_key))
rm(heavy_tips)
# Group members, kept for the links written after `resolved` is freed.
group_members <- resolved[!duplicated(resolved$sequence_id),
                          c("sequence_id", "group_id", "is_representative")]
rm(clono, clones, joined, resolved, db, heavy_seqs, frames, lineages, heavy_lineages, kept_rows)
invisible(gc())

# One lineage's work order; `slot` keeps output in build order, not finishing order.
slots <- 0L
unit_base <- function(key, build, part, label) {
  slots <<- slots + 1L
  gapped_tips <- tips_by_lineage[[key]]
  list(slot = slots, build = build, part = part, label = label,
       size = length(tip_rows[[key]]), lid = unname(real_lineage[key]),
       gapped_tips = if (is.null(gapped_tips)) character(0) else gapped_tips,
       gapped_germ = unname(gapped_germline[key]))
}
unit_for <- function(one, build, part, label) {
  c(unit_base(as.character(one$clone_id), build, part, label), list(one = one))
}
unit_rows <- function(key) {
  c(unit_base(key, "raxml", NULL, NULL), list(db = tips_db[tip_rows[[key]], , drop = FALSE]))
}

# IgPhyML fits one model per call, so each group is one call; if it fails, parsimony.
queue_igphyml <- function(subset, part, label) {
  if (!nrow(subset)) return(list(units = list(), finished = list()))
  progress(sprintf("Trees: IgPhyML on %d lineages", nrow(subset)))
  result <- attempt_build("igphyml", part, subset, fork_workers(threads))
  if (!usable_build(result) && !is.null(part)) {
    cat(sprintf("%s: partitioned build produced nothing (%s); rebuilding unpartitioned\n",
                label, why_build(result)))
    result <- attempt_build("igphyml", NULL, subset, fork_workers(threads))
  }
  if (usable_build(result)) {
    done <- lapply(seq_len(nrow(result)), function(i) {
      finish_lineage(list(tree = result[i, ], builder = "igphyml"),
                     unit_for(result[i, ], "igphyml", part, label))
    })
    trees_done <<- trees_done + sum(lengths(tip_rows[as.character(subset$clone_id)]))
    progress(sprintf("Trees: %.1f%%", 100 * trees_done / max(1, trees_total)))
    return(list(units = list(), finished = done))
  }
  cat(sprintf("%s: igphyml produced no tree (%s); falling back to maximum parsimony\n",
              label, why_build(result)))
  list(units = lapply(seq_len(nrow(subset)), function(i) unit_for(subset[i, ], "pratchet", NULL, label)),
       finished = list())
}

trees_total <- nrow(tips_db)
trees_started <- Sys.time()
progress("Trees: 0.0%")

formatted <- list(lineages = 0L, seqs = 0L, mixed = 0L, below_min = 0L)
queued <- list()
ig_keys <- lineage_keys[builder_of == "igphyml"]
if (length(ig_keys)) {
  ig_fmt <- tryCatch(format_lineages(tips_db[unlist(tip_rows[ig_keys], use.names = FALSE), , drop = FALSE],
                                     fork_workers(min(8L, fit_workers(threads, "formatClones")))),
                     error = function(e) { cat(sprintf("formatClones: %s\n", conditionMessage(e))); NULL })
  if (!is.null(ig_fmt) && nrow(ig_fmt)) {
    formatted$lineages <- nrow(ig_fmt)
    formatted$seqs <- sum(ig_fmt$seqs)
    if (min_tips > 2) {
      formatted$below_min <- sum(ig_fmt$seqs < min_tips)
      ig_fmt <- ig_fmt[ig_fmt$seqs >= min_tips, ]
    }
    ig_mixed <- if (chain == "HL") vapply(ig_fmt$data, function(x) any(x@locus != HEAVY), logical(1)) else rep(FALSE, nrow(ig_fmt))
    formatted$mixed <- sum(ig_mixed)
    queued[[1]] <- queue_igphyml(ig_fmt[ig_mixed, ], partition_for("igphyml"), "igphyml with light chains")
    queued[[2]] <- queue_igphyml(ig_fmt[!ig_mixed, ], NULL, "igphyml heavy chain only")
  }
  rm(ig_fmt)
}
units <- c(do.call(c, lapply(queued, function(q) q$units)),
           lapply(lineage_keys[builder_of == "raxml"], unit_rows))
finished <- do.call(c, lapply(queued, function(q) q$finished))
rm(queued, tips_db, tips_by_lineage)
invisible(gc())

if (length(units)) {
  tree_workers <- fit_pool(min(threads, length(units)))
  # Default outfile silences worker R output (captured anyway); child programs still print.
  cl <- parallel::makePSOCKcluster(tree_workers)
  parallel::clusterEvalQ(cl, {
    options(warn = 1)
    suppressMessages({library(dowser); library(ape)})
    NULL
  })
  parallel::clusterExport(cl, POOL_EXPORTS)
  finished <- c(finished, run_pool(cl, units))
  parallel::stopCluster(cl)
  rm(units)
}

# Merging and writing the results can take a while on large donors.
progress("Saving trees")
from_pool <- Filter(function(r) isTRUE(r$formatted), finished)
formatted$lineages <- formatted$lineages + length(from_pool)
formatted$seqs <- formatted$seqs + sum(vapply(from_pool, function(r) r$kept_seqs, integer(1)))
formatted$below_min <- formatted$below_min + sum(vapply(from_pool, function(r) r$below_min, logical(1)))
formatted$mixed <- formatted$mixed + sum(vapply(from_pool, function(r) r$mixed && !r$below_min, logical(1)))
stops_removed <- sum(vapply(finished, function(r) if (is.null(r$stops_removed)) 0L else r$stops_removed, integer(1)))
if (stops_removed > 0) {
  cat(sprintf("warning: %d sequence(s) with an inframe stop codon were removed. If you want to keep these sequences use the option filterstop=FALSE.\n",
              stops_removed))
}
cat(sprintf("formatClones(chain=%s) kept %d of %d distinct sequences (%.1f%%) in %d lineages\n",
            chain, formatted$seqs, trees_total, pct(formatted$seqs, trees_total), formatted$lineages))
if (min_tips > 2) {
  cat(sprintf("%d of %d lineages (%.1f%%) have fewer than %d tips and get no tree\n",
              formatted$below_min, formatted$lineages, pct(formatted$below_min, formatted$lineages), min_tips))
  if (formatted$below_min == formatted$lineages) finish_empty(sprintf("no lineage reached %d tips", min_tips), membership)
}
if (chain == "HL") {
  cat(sprintf("%d of %d lineages carry a light chain\n", formatted$mixed, formatted$lineages - formatted$below_min))
}

# Back in build order: IgPhyML then RAxML groups, light before heavy only, small lineages last.
GROUP_ORDER <- c("igphyml with light chains", "igphyml heavy chain only",
                 "raxml with light chains", "raxml heavy chain only")
finished <- Filter(function(r) !is.null(r$label), finished)
rank_of <- function(r) match(r$label, GROUP_ORDER, nomatch = length(GROUP_ORDER) + 1L)
finished <- finished[order(vapply(finished, rank_of, integer(1)),
                           vapply(finished, function(r) isTRUE(r$small), logical(1)),
                           vapply(finished, function(r) r$slot, integer(1)))]
for (label in intersect(GROUP_ORDER, vapply(finished, function(r) r$label, character(1)))) {
  own <- Filter(function(r) identical(r$label, label), finished)
  # The model each group asked for; a lineage that fell back says so in a note.
  part <- if (grepl("with light chains", label)) partition_for(sub(" .*", "", label)) else NULL
  cat(sprintf("%s: %d lineages, partition: %s\n", label, length(own), if (is.null(part)) "none" else part))
  small <- sum(vapply(own, function(r) isTRUE(r$small), logical(1)))
  if (small > 0) {
    cat(sprintf("%s: %d of %d lineages have fewer than 3 tips, below what FastTree and RAxML accept; building those with maximum parsimony\n",
                label, small, length(own)))
  }
  notes <- unlist(lapply(own, function(r) r$note))
  for (note in unique(notes)) cat(sprintf("%s: %d lineages %s\n", label, sum(notes == note), note))
}
built_ok <- Filter(function(r) !is.null(r$nodes), finished)
rm(finished)
if (!length(built_ok)) finish_empty("no lineage produced a tree", membership)
cat(sprintf("collapsed %d of %d internal nodes with identical reconstructed sequences\n",
            sum(vapply(built_ok, function(r) r$nodes_before - r$nodes_after, integer(1))),
            sum(vapply(built_ok, function(r) r$nodes_before, integer(1)))))
trees <- data.frame(clone_id = vapply(built_ok, function(r) r$clone_id, character(1)),
                    tree_builder = vapply(built_ok, function(r) r$builder, character(1)),
                    stringsAsFactors = FALSE)
anchor_rows <- Filter(Negate(is.null), lapply(built_ok, function(r) r$anchors))
nodes <- as.data.frame(dplyr::bind_rows(lapply(built_ok, function(r) r$nodes)))
rm(built_ok)
write.table(nodes, nodes_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

if (!is.null(anchor_path)) {
  distances <- if (length(anchor_rows)) do.call(rbind, anchor_rows) else
    as.data.frame(setNames(replicate(length(ANCHOR_COLUMNS), character(0), simplify = FALSE),
                           ANCHOR_COLUMNS))
  write.table(distances, anchor_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("anchor distances: %d clonotypes measured in %d anchored lineages\n",
              nrow(distances), length(anchor_rows)))
}

if (!is.null(builders_path)) {
  write.table(
    data.frame(
      lineage_id = unname(real_lineage[as.character(trees$clone_id)]),
      tree_builder = trees$tree_builder,
      stringsAsFactors = FALSE),
    builders_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

progress("Linking tips to clonotypes")
# Link each tip to every clonotype in its group; the representative is marked.
tips <- nodes[nodes$is_observed == "true" & nodes$label != "Germline",
              c("lineage_id", "node_id", "label")]
tip_group <- group_members$group_id[match(tips$label, group_members$sequence_id)]
links <- merge(
  data.frame(lineage_id = tips$lineage_id, node_id = tips$node_id, group_id = tip_group,
             stringsAsFactors = FALSE),
  group_members, by = "group_id")
links <- links[, c("lineage_id", "node_id", "sequence_id", "is_representative")]
links$is_representative <- ifelse(links$is_representative, "true", "false")
links <- links[order(links$lineage_id, links$node_id, links$sequence_id), , drop = FALSE]
links$link <- 1L
write.table(links, links_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("built %d trees, %d nodes, %d lineages skipped, tree building took %s\n",
            nrow(trees), nrow(nodes), length(skipped), clock(trees_started)))
progress("Trees: 100.0%")
close_log()
