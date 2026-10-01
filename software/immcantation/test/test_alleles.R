#!/usr/bin/env Rscript
# Checks for alleles.R. Run through test/run.sh.
# Inference may change the germline a sequence is compared to, never the sequence itself.
# Each run is the alignment step (trees.R --stage align) then alleles.R on its rows, as in the
# workflow. Inference is slow, so a few donors cover its routes; other runs stop before it.

work <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(work)) stop("usage: test_alleles.R <work dir>")
alleles_R <- "/app/alleles.R"
trees_R <- "/app/trees.R"
dir.create(work, recursive = TRUE, showWarnings = FALSE)

failures <- 0L; checks <- 0L
ok <- function(what, cond) {
  checks <<- checks + 1L
  if (isTRUE(cond)) cat(sprintf("  ok   %s\n", what))
  else { cat(sprintf("  FAIL %s\n", what)); failures <<- failures + 1L }
}
read_tsv <- function(p) read.delim(p, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
write_tsv <- function(d, p) write.table(d, p, sep = "\t", quote = FALSE, row.names = FALSE)
gene_of <- function(x) sub("[*].*$", "", sub(",.*$", "", x))
allele_of <- function(x) sub(",.*$", "", x)

# TIgGER's sample repertoire cut to a few alleles. It has no germline_alignment, so build it from the called allele.
db <- tigger::AIRRDb
gl <- tigger::SampleGermlineIGHV
names(gl) <- sub(",.*$", "", names(gl))
allele <- allele_of(db$v_call)
take <- c(
  # IGHV1-8*02 carries the novel IGHV1-8*02_G234T.
  "IGHV1-8*01" = 250, "IGHV1-8*02" = 250,
  # Two common alleles of one gene, for the mislabelled rows below.
  "IGHV1-2*02" = 250, "IGHV1-2*04" = 250,
  # Dropped by the genotype and must be put back: one too rare, one gene with no unmutated row.
  "IGHV1-2*03" = 7, "IGHV1-3*01" = 15,
  "IGHV1-18*01" = 100,
  # Blank germline, so absent from the reference; moves to other genes must be refused.
  "IGHV1-24*01" = 20)
rows <- unlist(lapply(names(take), function(a) head(which(allele == a), take[[a]])))
db <- db[rows, ]
allele <- allele[rows]
with_germline <- function(s, a) {
  n <- min(nchar(gl[[a]]), nchar(s))
  paste0(substr(gl[[a]], 1, n), substr(s, n + 1, nchar(s)))
}
# Three V mutations on each IGHV1-3 row, so none is unmutated.
mutated <- which(allele == "IGHV1-3*01")
for (i in mutated) for (at in c(60L, 120L, 180L)) {
  b <- substr(db$sequence_alignment[i], at, at)
  if (b %in% c("A", "C", "G", "T")) substr(db$sequence_alignment[i], at, at) <- if (b == "A") "C" else "A"
}
germ <- mapply(with_germline, db$sequence_alignment, allele, USE.NAMES = FALSE)
bait <- which(allele == "IGHV1-24*01")
germ[bait] <- strrep("-", nchar(germ[bait]))
donor <- data.frame(
  sequence_id = paste0("ck_", seq_len(nrow(db))),
  v_call = db$v_call, j_call = db$j_call, junction = db$junction,
  sequence_alignment = db$sequence_alignment, germline_alignment = germ,
  stringsAsFactors = FALSE)
# As merge leaves it: the dataset, and the assembled sequence the alignment join breaks ties on.
donor$dataset <- "0"
donor$main_sequence <- toupper(gsub("[.-]", "", donor$sequence_alignment))

# Call some *04 rows as *02, germline too; their bases still match *04, so inference should fix them.
onto <- "IGHV1-2*02"; from <- "IGHV1-2*04"
mislabelled <- head(which(allele == from), 100L)
donor$v_call[mislabelled] <- onto
donor$germline_alignment[mislabelled] <- vapply(donor$sequence_alignment[mislabelled],
  with_germline, character(1), a = onto, USE.NAMES = FALSE)

# V-side mismatches of rebuilt rows, counted apart from alleles.R: the first frame_left bases.
v_side_mut <- function(a) {
  vapply(seq_len(nrow(a)), function(i) {
    l <- as.integer(a$frame_left[i])
    if (is.na(l) || l == 0L) return(NA_integer_)
    q <- strsplit(toupper(substr(a$sequence_alignment[i], 1L, l)), "")[[1]]
    r <- strsplit(toupper(substr(a$germline_alignment[i], 1L, l)), "")[[1]]
    use <- q %in% c("A", "C", "G", "T") & r %in% c("A", "C", "G", "T")
    sum(q[use] != r[use])
  }, integer(1))
}

# The MiXCR route: a gene-level table without alignments, the alleles in an AIRR export.
airr_dir <- file.path(work, "airr-all")
dir.create(airr_dir, showWarnings = FALSE)
write_tsv(donor, file.path(airr_dir, "ds0__s1.tsv"))
mixcr_table <- donor
mixcr_table$v_call <- gene_of(mixcr_table$v_call)
mixcr_table$sequence_alignment <- ""
mixcr_table$germline_alignment <- ""
write_tsv(mixcr_table, file.path(work, "mixcr.tsv"))
write_tsv(donor, file.path(work, "donor.tsv"))
# As import-vdj-data hands it over: the gene in v_call, the allele in v_allele.
stripped <- donor
stripped$v_allele <- stripped$v_call
stripped$v_call <- gene_of(stripped$v_call)
write_tsv(stripped, file.path(work, "stripped.tsv"))
# Both in one donor, as merge leaves it: MiXCR rows get an empty v_allele and no alignment.
half <- which(seq_len(nrow(donor)) %% 2 == 0)
mixed <- stripped
mixed$v_allele[half] <- ""
mixed$sequence_alignment[half] <- ""
mixed$germline_alignment[half] <- ""
write_tsv(mixed, file.path(work, "mixed.tsv"))
# Half the rows lose 25 columns at the 5' end: still on one numbering, counted from the junction.
truncated <- donor
for (column in c("sequence_alignment", "germline_alignment")) {
  truncated[[column]][half] <- substr(truncated[[column]][half], 26, nchar(truncated[[column]][half]))
}
write_tsv(truncated, file.path(work, "truncated.tsv"))

# Inputs for the reference route, each stopping before inference.
calls <- c("sequence_id", "v_call", "j_call", "junction", "dataset")
write_tsv(donor[, calls], file.path(work, "bare.tsv"))
write_tsv(donor[seq_len(50), ], file.path(work, "shallow.tsv"))
# Half the germlines changed at six V positions, so rows disagree with their own allele.
disagreeing <- donor
for (i in half) {
  g <- strsplit(disagreeing$germline_alignment[i], "")[[1]]
  at <- which(g %in% c("A", "C", "G", "T"))[c(40, 80, 120, 160, 200, 240)]
  at <- at[!is.na(at)]
  g[at] <- ifelse(g[at] == "A", "C", "A")
  disagreeing$germline_alignment[i] <- paste(g, collapse = "")
}
write_tsv(disagreeing, file.path(work, "disagreeing.tsv"))
write_tsv(donor[0, ], file.path(work, "empty.tsv"))

run <- function(name, input, airr = NULL) {
  f <- function(x) file.path(work, paste0(name, x))
  extra <- if (is.null(airr)) character() else c("--airr-dir", airr)
  align_log <- system2("Rscript", c(trees_R, "--stage", "align", "--clonotypes", input,
                                    "--out-aligned", f("-aligned.tsv"), "--out-log", f("-align.log"),
                                    "--threads", "1", extra), stdout = TRUE, stderr = TRUE)
  log <- system2("Rscript", c(alleles_R, "--clonotypes", input, "--aligned", f("-aligned.tsv"),
                              "--align-log", f("-align.log"), "--out", f("-out.tsv"),
                              "--out-aligned", f("-aligned-out.tsv"), "--out-route", f("-route.json"),
                              "--threads", "1"), stdout = TRUE, stderr = TRUE)
  rd <- function(p) if (file.exists(p)) read_tsv(p)
  list(ok = is.null(attr(log, "status")) && is.null(attr(align_log, "status")),
       log = paste(c(align_log, log), collapse = "\n"),
       out = rd(f("-out.tsv")), before = rd(f("-aligned.tsv")), after = rd(f("-aligned-out.tsv")),
       route = if (file.exists(f("-route.json"))) jsonlite::fromJSON(f("-route.json")))
}
RUNS <- list(
  donor = list("donor", file.path(work, "donor.tsv")),
  stripped = list("stripped", file.path(work, "stripped.tsv")),
  mixcr = list("mixcr", file.path(work, "mixcr.tsv"), airr_dir),
  mixed = list("mixed", file.path(work, "mixed.tsv"), airr_dir),
  truncated = list("truncated", file.path(work, "truncated.tsv")),
  bare = list("bare", file.path(work, "bare.tsv")),
  shallow = list("shallow", file.path(work, "shallow.tsv")),
  disagreeing = list("disagreeing", file.path(work, "disagreeing.tsv")),
  empty = list("empty", file.path(work, "empty.tsv")))
results <- parallel::mclapply(RUNS, function(a) do.call(run, a),
                              mc.cores = length(RUNS), mc.preschedule = FALSE)
tail_of <- function(r) cat(substr(r$log, max(1, nchar(r$log) - 1500), nchar(r$log)), "\n")
heavy <- function(a) a[a$locus == "IGH", , drop = FALSE]
# Each donor row's call after the run, read off the rebuilt rows.
call_after <- function(r) heavy(r$after)$v_call[match(donor$sequence_id, heavy(r$after)$sequence_id)]

cat("== one donor through TIgGER, on the alignment step's rows ==\n")
r <- results$donor
ok("runs", isTRUE(r$ok))
if (isTRUE(r$ok) && !is.null(r$after)) {
  ok("the route taken is tigger", r$route$route == "tigger")
  ok("it pools the rebuilt rows and repeats what the alignment step did",
     grepl("from the aligned table", r$log, fixed = TRUE) &&
       grepl("alignment: IMGT gaps removed from", r$log, fixed = TRUE))
  b <- heavy(r$before); a <- heavy(r$after)
  ok("rows and their order survive", identical(a$sequence_id, b$sequence_id))
  ok("sequence_alignment, junction and j_call are never touched",
     identical(a$sequence_alignment, b$sequence_alignment) &&
       identical(a$junction, b$junction) && identical(a$j_call, b$j_call))
  ok("germline and query still line up", all(nchar(a$germline_alignment) == nchar(a$sequence_alignment)))
  ok("only the V side of a germline changes",
     all(substr(a$germline_alignment, as.integer(a$frame_left) + 1L, nchar(a$germline_alignment)) ==
           substr(b$germline_alignment, as.integer(b$frame_left) + 1L, nchar(b$germline_alignment))))
  now <- call_after(r)
  novel <- grepl("_", now, fixed = TRUE)
  kept <- regmatches(r$log, regexec("genotype of \\d+ allele\\(s\\) over \\d+ gene\\(s\\), (\\d+) of the", r$log))[[1]]
  ok("novel alleles reach the genotype", length(kept) == 2 && as.integer(kept[2]) > 0)
  ok("and are applied, each to sequences of its own gene",
     any(grepl("^IGHV1-8\\*02_", now)) && all(gene_of(now[novel]) == gene_of(donor$v_call[novel])))
  ok("the mislabelled sequences go back to the allele they came from",
     all(allele_of(now[mislabelled]) == from))
  ok("the clonotype table carries the same new calls",
     identical(r$out$v_call[match(a$sequence_id, r$out$sequence_id)][a$v_call != b$v_call],
               a$v_call[a$v_call != b$v_call]))
  before <- v_side_mut(b); after <- v_side_mut(a)
  rewritten <- which(a$germline_alignment != b$germline_alignment)
  ok("an offered allele no closer than the called one is refused",
     grepl("sits no closer to them than the one they were called against", r$log))
  ok("every rewritten row ends strictly closer to its germline",
     length(rewritten) > 0 && all(!is.na(after[rewritten]) & after[rewritten] < before[rewritten]))
  ok("and no row ends further from it", all(after <= before, na.rm = TRUE))
  ok("calls and germlines move together", setequal(rewritten, which(a$v_call != b$v_call)))
  ok("rows with no germline of their own keep their call", identical(now[bait], donor$v_call[bait]))
  ok("no sequence changes gene", all(gene_of(now) == gene_of(donor$v_call)))
  put_back <- regmatches(r$log, regexpr("put back: [^\n]*", r$log))
  ok("reference alleles the genotype dropped are put back",
     length(put_back) == 1 && grepl("IGHV1-2*03", put_back, fixed = TRUE) &&
       grepl("IGHV1-3*01", put_back, fixed = TRUE))
} else tail_of(r)

# The same donor reaching the rows another way: the call it groups by, and its novel allele.
for (name in c("stripped", "mixcr", "mixed", "truncated")) {
  cat(sprintf("== %s ==\n", name))
  r <- results[[name]]
  ok(sprintf("%s: TIgGER runs over the whole donor", name),
     isTRUE(r$ok) && identical(r$route$route, "tigger") &&
       grepl(sprintf("pooled %d heavy sequences", nrow(donor) - length(bait)), r$log, fixed = TRUE))
  if (!isTRUE(r$ok) || is.null(r$after)) { tail_of(r); next }
  now <- call_after(r)
  ok(sprintf("%s: the mislabelled sequences go back, and the novel allele is found", name),
     all(allele_of(now[mislabelled]) == from) && any(grepl("^IGHV1-8\\*02_", now)))
}

# Each case must keep the reference alleles, not fail, and pass both tables through.
for (case in list(
  list(name = "bare", says = "no clonotype has an alignment"),
  list(name = "shallow", says = "below the"),
  list(name = "disagreeing", says = "not on one numbering"),
  list(name = "empty", says = "no clonotypes for this donor"))) {
  cat(sprintf("== reference route: %s ==\n", case$name))
  r <- results[[case$name]]
  ok(sprintf("%s: runs, keeps the reference alleles, and says why", case$name),
     isTRUE(r$ok) && identical(r$route$route, "reference") && grepl(case$says, r$route$reason, fixed = TRUE))
  if (!isTRUE(r$ok)) { tail_of(r); next }
  original <- read_tsv(file.path(work, paste0(case$name, ".tsv")))
  ok(sprintf("%s: both tables pass through untouched", case$name),
     (identical(r$out, original) || (nrow(original) == 0 && nrow(r$out) == 0)) &&
       (identical(r$after, r$before) || (nrow(r$before) == 0 && nrow(r$after) == 0)))
}

cat(sprintf("\n%d checks, %d failures\n", checks, failures))
quit(status = if (failures > 0L) 1L else 0L)
