#!/usr/bin/env Rscript
# Checks for alleles.R. Run through test/run.sh.
# Inference may change the germline a sequence is compared to, never the sequence itself.
# Inference is slow, so one donor covers every route through it; other runs stop before it.

work <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(work)) stop("usage: test_alleles.R <work dir>")
alleles_R <- "/app/alleles.R"
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

# Call some *04 rows as *02, germline too; their bases still match *04, so inference should fix them.
onto <- "IGHV1-2*02"; from <- "IGHV1-2*04"
mislabelled <- head(which(allele == from), 100L)
donor$v_call[mislabelled] <- onto
donor$germline_alignment[mislabelled] <- vapply(donor$sequence_alignment[mislabelled],
  with_germline, character(1), a = onto, USE.NAMES = FALSE)

# V mismatches before the junction, counted apart from alleles.R so it checks it independently.
v_mut <- function(s, g, junction) {
  vapply(seq_along(s), function(i) {
    sc <- strsplit(s[i], "")[[1]]
    gc <- strsplit(g[i], "")[[1]]
    if (length(sc) != length(gc)) return(NA_integer_)
    bases <- which(!(sc %in% c(".", "-")))
    at <- regexpr(toupper(junction[i]), toupper(paste(sc[bases], collapse = "")), fixed = TRUE)
    if (at < 0 || bases[at] <= 1) return(NA_integer_)
    cols <- seq_len(bases[at] - 1L)
    q <- toupper(sc[cols]); r <- toupper(gc[cols])
    use <- q %in% c("A", "C", "G", "T") & r %in% c("A", "C", "G", "T")
    sum(q[use] != r[use])
  }, integer(1))
}

# Half the rows keep alignments on the table, half in an AIRR file; some table rows share a join key with AIRR rows.
join_key <- paste(gene_of(donor$v_call), gene_of(donor$j_call), donor$junction)
in_airr <- which(seq_len(nrow(donor)) %% 2 == 0)
in_table <- setdiff(seq_len(nrow(donor)), in_airr)
table_in <- donor
table_in$sequence_alignment[in_airr] <- ""
table_in$germline_alignment[in_airr] <- ""
airr_dir <- file.path(work, "airr-in")
dir.create(airr_dir, showWarnings = FALSE)
write_tsv(donor[in_airr, ], file.path(airr_dir, "ds0__s1.tsv"))
write_tsv(table_in, file.path(work, "donor.tsv"))

# Inputs for the reference route, each stopping before inference.
calls <- c("sequence_id", "v_call", "j_call", "junction")
write_tsv(donor[, calls], file.path(work, "bare.tsv"))
bare_airr <- file.path(work, "bare-airr")
dir.create(bare_airr, showWarnings = FALSE)
write_tsv(donor[in_airr, calls], file.path(bare_airr, "ds0__s1.tsv"))
write_tsv(donor[seq_len(50), ], file.path(work, "shallow.tsv"))
# Every other row starts 25 columns in, so rows no longer share one numbering.
misnumbered <- donor
for (column in c("sequence_alignment", "germline_alignment")) {
  misnumbered[[column]][in_airr] <- substr(misnumbered[[column]][in_airr], 26,
                                           nchar(misnumbered[[column]][in_airr]))
}
write_tsv(misnumbered, file.path(work, "misnumbered.tsv"))
write_tsv(donor[0, ], file.path(work, "empty.tsv"))
# As import-vdj-data hands it over: the gene in v_call, the allele in v_allele.
stripped <- donor
stripped$v_allele <- stripped$v_call
stripped$v_call <- gene_of(stripped$v_call)
write_tsv(stripped, file.path(work, "stripped.tsv"))
# MiXCR: a gene-level table without alignments, the alleles in the AIRR exports.
# The whole donor, as half of it is too few for TIgGER.
all_airr <- file.path(work, "airr-all")
dir.create(all_airr, showWarnings = FALSE)
write_tsv(donor, file.path(all_airr, "ds0__s1.tsv"))
mixcr_table <- donor
mixcr_table$v_call <- gene_of(mixcr_table$v_call)
mixcr_table$sequence_alignment <- ""
mixcr_table$germline_alignment <- ""
write_tsv(mixcr_table, file.path(work, "mixcr.tsv"))
# Both in one donor, as merge leaves it: MiXCR rows get an empty v_allele.
mixed <- stripped
mixed$v_allele[in_airr] <- ""
mixed$sequence_alignment[in_airr] <- ""
mixed$germline_alignment[in_airr] <- ""
write_tsv(mixed, file.path(work, "mixed.tsv"))

run <- function(name, input, airr_dir = NULL) {
  out <- file.path(work, paste0(name, "-out.tsv"))
  route <- file.path(work, paste0(name, "-route.json"))
  out_airr <- file.path(work, paste0(name, "-airr-out"))
  extra <- if (is.null(airr_dir)) character() else c("--airr-dir", airr_dir, "--out-airr-dir", out_airr)
  log <- system2("Rscript", c(alleles_R, "--clonotypes", input, "--out", out,
                              "--out-route", route, "--threads", "1", extra),
                 stdout = TRUE, stderr = TRUE)
  airr_out <- if (!is.null(airr_dir) && dir.exists(out_airr)) {
    setNames(lapply(sort(list.files(out_airr, full.names = TRUE)), read_tsv), sort(list.files(out_airr)))
  }
  list(ok = is.null(attr(log, "status")), log = paste(log, collapse = "\n"),
       out = if (file.exists(out)) read_tsv(out), airr = airr_out,
       route = if (file.exists(route)) jsonlite::fromJSON(route))
}
RUNS <- list(
  donor = list("donor", file.path(work, "donor.tsv"), airr_dir),
  bare = list("bare", file.path(work, "bare.tsv"), bare_airr),
  shallow = list("shallow", file.path(work, "shallow.tsv")),
  misnumbered = list("misnumbered", file.path(work, "misnumbered.tsv")),
  empty = list("empty", file.path(work, "empty.tsv")),
  stripped = list("stripped", file.path(work, "stripped.tsv")),
  mixcr = list("mixcr", file.path(work, "mixcr.tsv"), all_airr),
  mixed = list("mixed", file.path(work, "mixed.tsv"), airr_dir))
results <- parallel::mclapply(RUNS, function(a) do.call(run, a),
                              mc.cores = length(RUNS), mc.preschedule = FALSE)
tail_of <- function(r) cat(substr(r$log, max(1, nchar(r$log) - 1500), nchar(r$log)), "\n")

cat("== one donor through TIgGER: table rows and AIRR rows ==\n")
r <- results$donor
ok("runs", isTRUE(r$ok))
back <- r$airr[["ds0__s1.tsv"]]
if (isTRUE(r$ok) && !is.null(back)) {
  ok("the route taken is tigger", r$route$route == "tigger")
  ok("the reference comes from the data, pooled from the table and the AIRR file",
     grepl("germline database derived from the data", r$log) &&
       grepl("from the clonotype table and 1 AIRR file", r$log))
  # Each row's current state, from the table or the AIRR file.
  airr_row <- match(donor$sequence_id, back$sequence_id)
  now <- data.frame(v_call = ifelse(is.na(airr_row), r$out$v_call, back$v_call[airr_row]),
                    sequence_alignment = ifelse(is.na(airr_row), r$out$sequence_alignment,
                                                back$sequence_alignment[airr_row]),
                    germline_alignment = ifelse(is.na(airr_row), r$out$germline_alignment,
                                                back$germline_alignment[airr_row]),
                    stringsAsFactors = FALSE)
  ok("rows and their order survive, in the table and in the AIRR file",
     identical(r$out$sequence_id, donor$sequence_id) &&
       identical(back$sequence_id, donor$sequence_id[in_airr]))
  ok("sequence_alignment, junction and j_call are never touched",
     identical(now$sequence_alignment, donor$sequence_alignment) &&
       identical(r$out$junction, donor$junction) && identical(r$out$j_call, donor$j_call))
  ok("germline and query still line up",
     all(nchar(now$germline_alignment) == nchar(now$sequence_alignment)))

  novel <- grepl("_", now$v_call, fixed = TRUE)
  kept <- regmatches(r$log, regexec("genotype of \\d+ allele\\(s\\) over \\d+ gene\\(s\\), (\\d+) of the", r$log))[[1]]
  ok("novel alleles reach the genotype", length(kept) == 2 && as.integer(kept[2]) > 0)
  ok("and are applied, each to sequences of its own gene",
     "IGHV1-8*02_G234T" %in% now$v_call &&
       all(gene_of(now$v_call[novel]) == gene_of(donor$v_call[novel])))
  ok("the mislabelled sequences go back to the allele they came from",
     all(allele_of(now$v_call[mislabelled]) == from))
  ok("in the table and in the AIRR file alike",
     any(mislabelled %in% in_table) && any(mislabelled %in% in_airr))

  # A move must bring the row strictly closer and stay in its gene.
  before <- v_mut(donor$sequence_alignment, donor$germline_alignment, donor$junction)
  after <- v_mut(now$sequence_alignment, now$germline_alignment, donor$junction)
  rewritten <- which(now$germline_alignment != donor$germline_alignment)
  ok("an offered allele no closer than the called one is refused",
     grepl("sits no closer to them than the one they were called against", r$log))
  ok("every rewritten row ends strictly closer to its germline",
     length(rewritten) > 0 && all(!is.na(after[rewritten]) & after[rewritten] < before[rewritten]))
  ok("and no row ends further from it", all(after <= before, na.rm = TRUE))
  ok("calls and germlines move together",
     setequal(rewritten, which(now$v_call != donor$v_call)))
  shared <- in_table[join_key[in_table] %in% join_key[in_airr]]
  ok("a table row sharing a join key with an AIRR row keeps its own call",
     length(shared) > 0 && identical(r$out$v_call[shared] != donor$v_call[shared],
                                     r$out$germline_alignment[shared] != donor$germline_alignment[shared]))
  moved <- regmatches(r$log, regexec("(\\d+) sequences keep their original call: the reassignment would have moved them to another gene", r$log))[[1]]
  ok("a move to another gene is refused, for exactly the rows with no reference of their own",
     length(moved) == 2 && as.integer(moved[2]) == length(bait))
  ok("no sequence changes gene", all(gene_of(now$v_call) == gene_of(donor$v_call)))
  # Dropped alleles are put back so their rows are not pushed onto another gene.
  put_back <- regmatches(r$log, regexpr("put back: [^\n]*", r$log))
  ok("reference alleles the genotype dropped are put back",
     length(put_back) == 1 && grepl("IGHV1-2*03", put_back, fixed = TRUE) &&
       grepl("IGHV1-3*01", put_back, fixed = TRUE))

  # The align stage joins on V gene, J gene and junction, so both sides must move together.
  # Only unique keys are checked; shared keys are ambiguous by design.
  unique_key <- !(duplicated(join_key[in_airr]) | duplicated(join_key[in_airr], fromLast = TRUE))
  ok("the clonotype table takes its AIRR rows' new calls, so the join still matches",
     any(back$v_call != donor$v_call[in_airr]) &&
       identical(r$out$v_call[in_airr][unique_key], back$v_call[unique_key]))
} else tail_of(r)

cat("== gene-level v_call with the allele beside it ==\n")
r <- results$stripped
ok("stripped: groups by v_allele, so the alignments are on one numbering and TIgGER runs",
   isTRUE(r$ok) && r$route$route == "tigger")
if (!isTRUE(r$ok) || r$route$route != "tigger") tail_of(r)

# MiXCR rows reach TIgGER through the AIRR exports; their gene-level table rows follow by join key.
for (name in c("mixcr", "mixed")) {
  cat(sprintf("== %s ==\n", name))
  r <- results[[name]]
  back <- r$airr[["ds0__s1.tsv"]]
  airr_rows <- if (name == "mixcr") seq_len(nrow(donor)) else in_airr
  ok(sprintf("%s: TIgGER runs, pooled from the table and the AIRR file", name),
     isTRUE(r$ok) && r$route$route == "tigger" &&
       grepl(sprintf("pooled %d heavy sequences with alignments from the clonotype table and 1 AIRR file",
                     nrow(donor)), r$log, fixed = TRUE))
  if (!isTRUE(r$ok) || is.null(back)) { tail_of(r); next }
  fixed_airr <- intersect(mislabelled, airr_rows)
  ok(sprintf("%s: mislabelled AIRR rows go back to their allele", name),
     length(fixed_airr) > 0 && all(allele_of(back$v_call[match(fixed_airr, airr_rows)]) == from))
  unique_key <- !(duplicated(join_key[airr_rows]) | duplicated(join_key[airr_rows], fromLast = TRUE))
  # Rows whose AIRR call did not change keep the gene-level call they came with.
  # Unique join keys only: a shared key takes whichever of its rows moved.
  changed <- back$v_call != donor$v_call[airr_rows]
  input <- read_tsv(file.path(work, paste0(name, ".tsv")))
  ok(sprintf("%s: MiXCR table rows take their AIRR row's new call, the rest keep theirs", name),
     any(changed & unique_key) &&
       identical(r$out$v_call[airr_rows][changed & unique_key], back$v_call[changed & unique_key]) &&
       identical(r$out$v_call[airr_rows][!changed & unique_key], input$v_call[airr_rows][!changed & unique_key]))
  if (name == "mixed") {
    fixed_table <- intersect(mislabelled, in_table)
    ok("mixed: mislabelled imported rows go back to their allele",
       length(fixed_table) > 0 && all(allele_of(r$out$v_call[fixed_table]) == from))
  }
}

# Each case must keep the reference alleles, not fail.
for (case in list(
  list(name = "bare", says = "no table carries", input = "bare.tsv"),
  list(name = "shallow", says = "below the", input = "shallow.tsv"),
  list(name = "misnumbered", says = "not on one numbering", input = "misnumbered.tsv"),
  list(name = "empty", says = "no clonotypes for this donor", input = "empty.tsv"))) {
  cat(sprintf("== reference route: %s ==\n", case$name))
  r <- results[[case$name]]
  ok(sprintf("%s: runs, keeps the reference alleles, and says why", case$name),
     isTRUE(r$ok) && r$route$route == "reference" && grepl(case$says, r$route$reason, fixed = TRUE))
  if (!isTRUE(r$ok)) { tail_of(r); next }
  original <- read_tsv(file.path(work, case$input))
  ok(sprintf("%s: the table passes through untouched", case$name),
     identical(r$out, original) || (nrow(original) == 0 && nrow(r$out) == 0))
  if (case$name == "bare") {
    ok("bare: the AIRR files pass through untouched",
       identical(r$airr[["ds0__s1.tsv"]], read_tsv(file.path(bare_airr, "ds0__s1.tsv"))))
  }
}

cat(sprintf("\n%d checks, %d failures\n", checks, failures))
quit(status = if (failures > 0L) 1L else 0L)
