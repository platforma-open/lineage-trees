#!/usr/bin/env Rscript
# Infers a donor's V alleles with TIgGER, pooled over all the donor's datasets, against
# germlines taken from the input's own alignments (grouped by v_call), not a packaged set.
# Reads the alignment step's rows, rebuilt around the junction, so positions agree whatever
# tool or gapping produced them. Rewrites heavy-chain v_call and the V part of
# germline_alignment in those rows; never the sequence.
# Any failure leaves the donor on the reference alleles and writes the reason to --out-route.

suppressMessages({library(tigger)})

# Helpers shared with the other stage script, next to this one.
local({
  me <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(normalizePath(me)), "common.R"))
})

clonotypes_path <- opt("--clonotypes")
out_path <- opt("--out")
# The alignment step's rows: each rebuilt around its junction, V side first (`frame_left` bases).
aligned_path <- opt("--aligned")
out_aligned_path <- opt("--out-aligned")
align_log_path <- opt("--align-log", required = FALSE)
route_path <- opt("--out-route", required = FALSE)
threads <- {
  t <- opt("--threads", required = FALSE)
  if (is.null(t)) 1L else max(1L, as.integer(t))
}
# Minimum pooled sequences per donor. TIgGER finds novel alleles from how mutation
# at each position rises with overall load, which needs many sequences.
min_sequences <- {
  m <- opt("--min-sequences", required = FALSE)
  if (is.null(m)) 200L else max(1L, as.integer(m))
}

HEAVY <- "IGH"
# The last V positions are shaped by trimming, so polymorphisms there are mostly artefacts.
THREE_PRIME_MARGIN <- 10L
# TIgGER takes positions 1 to 312 as the V (IMGT's length), whatever its input.
IMGT_V_LENGTH <- 312L
# Sequences per allele given to the novel allele search.
MAX_PER_ALLELE <- 50000L
ALIGNED_NEEDED <- c("sequence_id", "v_call", "j_call", "junction", "sequence_alignment",
                    "germline_alignment", "locus", "frame_left")

# Flushes, so a killed run's log ends at the step that was running.
say <- function(...) { cat(sprintf(...)); flush(stdout()) }

write_tsv <- function(d, p) write.table(d, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

allele_of <- function(x) sub(",.*$", "", x)
present <- function(x) !is.na(x) & nzchar(x)

progress("Reading alignments")
clono <- read_tsv(clonotypes_path)
aligned <- if (file.exists(aligned_path)) read_tsv(aligned_path) else data.frame()

# Always writes both tables; this stage never fails the run.
finish <- function(route, reason) {
  cat(sprintf("%s\n", reason))
  write_tsv(clono, out_path)
  write_tsv(aligned, out_aligned_path)
  if (!is.null(route_path)) {
    # TIgGER's error messages span lines, so the reason is JSON-escaped.
    writeLines(jsonlite::toJSON(list(route = route, reason = reason), auto_unbox = TRUE), route_path)
  }
  quit(status = 0)
}

# What the alignment step did to the rows, shown here because its own log is not.
if (!is.null(align_log_path) && file.exists(align_log_path)) {
  notes <- grep("IMGT gaps removed|lost their alignment|insertion columns dropped",
                readLines(align_log_path), value = TRUE)
  for (n in notes) cat(sprintf("alignment: %s\n", trimws(sub("^\\[[^]]*\\]\\s*", "", n))))
}

if (!nrow(clono)) finish("reference", "no clonotypes for this donor")
if (!nrow(aligned) || !all(ALIGNED_NEEDED %in% names(aligned))) {
  finish("reference", "no clonotype has an alignment, so there is nothing to infer from")
}

# The V call to group by: an import's `v_allele` where it carries one beside a gene-level
# `v_call`, else the aligned row's call (MiXCR's AIRR rows carry the allele).
call_of <- aligned$v_call
if ("v_allele" %in% names(clono)) {
  own <- clono$v_allele[match(aligned$sequence_id, clono$sequence_id)]
  call_of <- ifelse(present(own), own, call_of)
}

progress("Pooling heavy sequences")
# The pool: every heavy row with a V side. Positions count from the junction, so they mean
# the same germline position whatever tool, gapping or start the alignment had.
left <- suppressWarnings(as.integer(aligned$frame_left))
rows <- which(aligned$locus == HEAVY & present(aligned$sequence_alignment) &
              present(aligned$germline_alignment) & present(call_of) &
              !is.na(left) & left > 0L & substr(gene_of(call_of), 1, 3) == HEAVY)
pool <- data.frame(
  row = rows,
  id = paste0("r", rows),
  sequence_id = aligned$sequence_id[rows],
  v_call = call_of[rows],
  j_call = aligned$j_call[rows],
  junction = aligned$junction[rows],
  v_query = toupper(substr(aligned$sequence_alignment[rows], 1L, left[rows])),
  v_ref = toupper(substr(aligned$germline_alignment[rows], 1L, left[rows])),
  stringsAsFactors = FALSE)
pool$allele <- allele_of(pool$v_call)
# One vote per distinct rearrangement: one heavy chain read as several clonotypes (paired with
# two light chains, or the same clone in two datasets) would otherwise count each time. J and
# junction stay in the key, since TIgGER reads their variety as evidence for a novel allele.
pool$key <- paste(pool$allele, pool$j_call, pool$junction, pool$v_query, sep = "\r")
distinct <- !duplicated(pool$key)
cat(sprintf("pooled %d heavy sequences with a V side from the aligned table, %d of them distinct\n",
            nrow(pool), sum(distinct)))
if (sum(distinct) < min_sequences) {
  finish("reference", sprintf("%d distinct heavy sequences with alignments, below the %d TIgGER needs to tell a novel allele from noise",
                              sum(distinct), min_sequences))
}

# Too short a stretch of germline to tell one allele of a gene from another.
readable <- nchar(pool$v_ref) >= 30L
if (sum(readable) < min_sequences) {
  finish("reference", sprintf("only %d of %d sequences have 30 or more germline bases before the junction, below the %d needed",
                              sum(readable), nrow(pool), min_sequences))
}

# Left-padded to one width, so position i is the same distance from the junction in every row.
width <- max(nchar(pool$v_ref[readable]))
pad_left <- function(x) paste0(strrep(".", width - nchar(x)), x)

# The germline database, derived: a row's V germline is the allele it was called against.
# Per allele the longest is kept; shorter reads cover less of its 5' end.
by_allele <- split(pool$v_ref[readable], pool$allele[readable])
germline_db <- vapply(by_allele, function(v) v[which.max(nchar(v))], character(1))
# Rows that disagree with their allele where both have bases are not on one numbering.
suffix_agrees <- function(a, b) {
  n <- min(nchar(a), nchar(b))
  n >= 30 && identical(substr(a, nchar(a) - n + 1L, nchar(a)), substr(b, nchar(b) - n + 1L, nchar(b)))
}
conflicts <- vapply(names(by_allele), function(a) {
  sum(!vapply(by_allele[[a]], suffix_agrees, logical(1), germline_db[[a]]))
}, integer(1))
span <- nchar(germline_db)
cat(sprintf("V germline bases before the junction, per allele: median %d, range %d to %d\n",
            as.integer(median(span)), min(span), max(span)))
cat(sprintf("germline database derived from the data: %d allele(s) over %d gene(s), %d of %d sequences disagree with their own allele\n",
            length(germline_db), length(unique(gene_of(names(germline_db)))),
            sum(conflicts), sum(readable)))
if (sum(conflicts) > 0.1 * sum(readable)) {
  finish("reference", sprintf("%d of %d sequences disagree with the germline of the allele they were called against, so the alignments are not on one numbering",
                              sum(conflicts), sum(readable)))
}
if (length(germline_db) < 2) {
  finish("reference", "fewer than two V alleles were called, so there is nothing to tell apart")
}

# TIgGER reads "." as an empty position, so a germline short at its 5' end is padded there.
# The scan stops short of the junction by the 3' margin.
short <- sum(nchar(germline_db) < width)
germline_db <- vapply(germline_db, pad_left, character(1))
pos_range <- seq_len(width - THREE_PRIME_MARGIN)
cat(sprintf("germlines aligned at the junction over %d positions (%d of %d shorter at the 5' end); scanning positions 1 to %d\n",
            width, short, length(germline_db), max(pos_range)))

# Genotyping needs genes with several called alleles. A `*00` library gives one per
# gene; the run still goes on, since novel allele search does not need a choice.
per_gene <- table(gene_of(names(germline_db)))
unresolved <- sum(grepl("\\*00$", names(germline_db)))
if (max(per_gene) == 1) {
  cat(sprintf("every gene has exactly one called allele (%d of %d named *00), so the genotype step has nothing to choose between; only novel alleles can be found\n",
              unresolved, length(germline_db)))
} else {
  cat(sprintf("alleles per gene: median %g, max %d (%d of %d named *00)\n",
              median(per_gene), max(per_gene), unresolved, length(germline_db)))
}

# The V alone, aligned at the junction and padded past its end: TIgGER's 312 cut then
# ends in padding every row shares rather than in the CDR3.
once <- pool[distinct, , drop = FALSE]
db <- data.frame(
  sequence_id = once$id,
  v_call = once$v_call,
  j_call = once$j_call,
  junction = once$junction,
  junction_length = nchar(once$junction),
  sequence_alignment = paste0(pad_left(once$v_query), strrep(".", max(0L, IMGT_V_LENGTH - width))),
  stringsAsFactors = FALSE)

# Input shape, logged to help debug runs killed inside `findNovelAlleles`.
per_allele <- sort(table(once$allele), decreasing = TRUE)
say("per allele: %d alleles, largest %d sequences, median %g, smallest %d\n",
    length(per_allele), max(per_allele), median(per_allele), min(per_allele))
say("alignment width: median %d, range %d to %d; junction length median %d, range %d to %d\n",
    as.integer(median(nchar(db$sequence_alignment))), min(nchar(db$sequence_alignment)),
    max(nchar(db$sequence_alignment)), as.integer(median(db$junction_length)),
    min(db$junction_length), max(db$junction_length))

# TIgGER's memory follows the largest allele's sequence count: one allele of 435,839 used
# 14.5 GiB. A seeded sample per allele bounds it and keeps far more than its 200-sequence floors.
set.seed(1L)
search_rows <- unlist(lapply(split(seq_len(nrow(db)), once$allele), function(i) {
  if (length(i) > MAX_PER_ALLELE) sample(i, MAX_PER_ALLELE) else i
}), use.names = FALSE)
capped <- sum(per_allele > MAX_PER_ALLELE)
if (capped) say("%d allele(s) sampled down to %d sequences for the novel allele search\n",
                capped, MAX_PER_ALLELE)

progress("Looking for novel alleles")
say("looking for novel alleles over %d sequences, single process\n", length(search_rows))
novel <- tryCatch(findNovelAlleles(db[sort(search_rows), , drop = FALSE], germline_db,
                                   pos_range = pos_range, nproc = 1L),
                  error = function(e) e)
if (inherits(novel, "error")) {
  finish("reference", sprintf("TIgGER could not look for novel alleles: %s", conditionMessage(novel)))
}
novel_rows <- selectNovel(novel)
cat(sprintf("TIgGER: %d novel allele(s) found\n", nrow(novel_rows)))

# With one allele per gene and no novel allele, reassignment can only move a sequence to
# another gene, which corrupts its germline and clustering class. Checked after the search.
with_novel <- table(gene_of(unique(c(names(germline_db), as.character(novel_rows$polymorphism_call)))))
if (max(with_novel) == 1) {
  finish("reference", sprintf("every gene has one allele and no novel allele was found, so a reassignment could only move a sequence to a different gene rather than refine it within its own; %d allele(s) over %d gene(s)",
                              length(germline_db), length(per_gene)))
}

progress("Inferring the genotype")
say("inferring the genotype\n")
# `find_unmutated` must be on: only then does `inferGenotype` add the novel alleles. It
# can also drop reference alleles, which are put back below. With no unmutated sequences
# it stops, and the donor keeps its reference alleles.
genotype <- tryCatch(inferGenotype(db, germline_db = germline_db, novel = novel, find_unmutated = TRUE),
                     error = function(e) e)
if (inherits(genotype, "error")) {
  finish("reference", sprintf("TIgGER could not infer a genotype: %s", conditionMessage(genotype)))
}
genotype_db <- tryCatch(genotypeFasta(genotype, germline_db, novel),
                        error = function(e) e)
if (inherits(genotype_db, "error") || !length(genotype_db)) {
  finish("reference", "TIgGER's genotype held no alleles")
}
novel_found <- as.character(novel_rows$polymorphism_call)
novel_kept <- intersect(novel_found, names(genotype_db))
cat(sprintf("TIgGER: genotype of %d allele(s) over %d gene(s), %d of the %d novel allele(s) found among them\n",
            length(genotype_db), length(unique(gene_of(names(genotype_db)))),
            length(novel_kept), length(novel_found)))
# Per gene: kept alleles with counts, and novel alleles left out.
if (all(c("gene", "alleles", "counts") %in% names(genotype))) {
  for (i in seq_len(nrow(genotype))) {
    g <- genotype$gene[i]
    kept_alleles <- strsplit(genotype$alleles[i], ",")[[1]]
    counts <- strsplit(as.character(genotype$counts[i]), ",")[[1]]
    dropped <- setdiff(novel_found[gene_of(novel_found) == g], paste0(g, "*", kept_alleles))
    say("  %s: %s%s\n", g,
        paste(sprintf("%s (%s)", kept_alleles, counts[seq_along(kept_alleles)]), collapse = ", "),
        if (length(dropped)) paste0("; novel left out: ", paste(sub("^.*\\*", "", dropped), collapse = ", ")) else "")
  }
}
# Put back reference alleles the genotype dropped: `find_unmutated` drops them for being
# mutated, not absent. Novel alleles add to the reference, never replace it.
put_back <- setdiff(names(germline_db), names(genotype_db))
if (length(put_back)) {
  genotype_db <- c(genotype_db, germline_db[put_back])
  say("%d reference allele(s) put back that the genotype dropped, so every reference allele stays available\n",
      length(put_back))
  say("  put back: %s\n", paste(put_back, collapse = ", "))
}

progress("Reassigning alleles")
say("reassigning alleles\n")
reassigned <- tryCatch(reassignAlleles(db, genotype_db), error = function(e) e)
if (inherits(reassigned, "error")) {
  finish("reference", sprintf("TIgGER could not reassign alleles: %s", conditionMessage(reassigned)))
}
# The result is in `v_call_genotyped`; `v_call` is left unchanged.
if (!"v_call_genotyped" %in% names(reassigned)) {
  finish("reference", "TIgGER returned no v_call_genotyped, so nothing was reassigned")
}
# Each row takes the call of the distinct sequence it shares.
voter <- once$id[match(pool$key, once$key)]
new_call <- reassigned$v_call_genotyped[match(voter, reassigned$sequence_id)]
kept <- present(new_call) & !is.na(new_call)
# Offered: rows whose new allele differs from the original one.
offered <- kept & allele_of(new_call) != allele_of(pool$v_call)
cat(sprintf("TIgGER: %d of %d heavy sequences reassigned, %d offered a different allele, %d given back the allele they had\n",
            sum(kept), nrow(pool), sum(offered), sum(kept & !offered)))

# The germline follows the call over the row's own V side, counted from the junction.
# Where the allele is not known that far 5', the row keeps its own base. Junction and J stay.
follow_allele <- function(own, allele) {
  a <- strsplit(substr(allele, width - nchar(own) + 1L, width), "")[[1]]
  o <- strsplit(own, "")[[1]]
  a[a == "."] <- o[a == "."]
  paste(a, collapse = "")
}

# Mismatches per row, counting only columns where both sides have A, C, G or T.
# Done in blocks over a byte matrix for speed.
BASES <- utf8ToInt("ACGT")
is_base <- function(x) x == BASES[1L] | x == BASES[2L] | x == BASES[3L] | x == BASES[4L]
mismatches <- function(q, g) {
  out <- rep(NA_integer_, length(q))
  at <- which(!is.na(q) & !is.na(g) & nzchar(q))
  block <- 4096L
  for (from in seq_len(ceiling(length(at) / block))) {
    rows <- at[seq((from - 1L) * block + 1L, min(from * block, length(at)))]
    # Padded with spaces, which never count.
    width <- max(nchar(q[rows]), nchar(g[rows]))
    flatten <- function(x) as.integer(charToRaw(
      paste0(toupper(x), strrep(" ", width - nchar(x)), collapse = "")))
    a <- flatten(q[rows])
    b <- flatten(g[rows])
    out[rows] <- as.integer(colSums(matrix(is_base(a) & is_base(b) & a != b, nrow = width)))
  }
  out
}

progress("Applying reassignments")
# Which reassignments to accept. A move to another gene is refused: it would inflate
# mutations and change the clustering class, and the aligner chose the gene.
usable_call <- offered & !vapply(genotype_db[allele_of(new_call)], is.null, logical(1))
cross_gene <- usable_call & gene_of(new_call) != gene_of(pool$v_call)

# A reassignment must also bring the sequence strictly closer to its germline; a tie
# keeps the original call. Needed because the genotype can leave a gene with only
# distant novel alleles.
proposed <- rep(NA_character_, nrow(pool))
comparable <- which(usable_call)
proposed[comparable] <- mapply(follow_allele, pool$v_ref[comparable],
                               unname(genotype_db[allele_of(new_call[comparable])]), USE.NAMES = FALSE)
mut_before <- mismatches(pool$v_query, pool$v_ref)
mut_after <- mismatches(pool$v_query, proposed)
not_closer <- usable_call & !cross_gene & !is.na(mut_before) & !is.na(mut_after) &
              mut_after >= mut_before

take <- which(usable_call & !cross_gene & !not_closer)
changed_gene <- sum(gene_of(new_call[take]) != gene_of(pool$v_call[take]))
if (any(cross_gene)) {
  say("%d sequences keep their original call: the reassignment would have moved them to another gene\n",
      sum(cross_gene))
}
if (any(not_closer)) {
  say("%d of the offered sequences keep their original call: the allele the genotype offered sits no closer to them than the one they were called against, so the reassignment would have added mutations rather than removed them (%d V mutations against %d)\n",
      sum(not_closer), sum(mut_after[not_closer]), sum(mut_before[not_closer]))
}
# The V side of each taken row's germline is replaced; the rest of the row is untouched.
rows_taken <- pool$row[take]
before <- aligned$germline_alignment[rows_taken]
aligned$germline_alignment[rows_taken] <- paste0(
  proposed[take], substr(before, left[rows_taken] + 1L, nchar(before)))
aligned$v_call[rows_taken] <- new_call[take]
rewritten <- length(take)
moved <- sum(before != aligned$germline_alignment[rows_taken])
# The clonotype table names the call too, for the tree step and the node tables.
at <- match(pool$sequence_id[take], clono$sequence_id)
clono$v_call[at[!is.na(at)]] <- new_call[take][!is.na(at)]

cat(sprintf("germlines rewritten on the V side for %d rows, %d of them changed, %d moved gene\n",
            rewritten, moved, changed_gene))

finish("tigger", sprintf("TIgGER inferred %d allele(s) for this donor, %d of the %d novel allele(s) found among them; %d sequence(s) took a new allele",
                         length(genotype_db), length(novel_kept), length(novel_found), length(take)))
