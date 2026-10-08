# R port of the bash raxml-ng-fasttree wrapper, for Windows (no bash). trees.R runs it through
# raxml-ng-fasttree.cmd there. Keep the two in step.
# No mode flag: FastTree draft, then raxml-ng --evaluate on it. Any mode flag: passed through.
# Programs come from RAXML_NG_EXEC and FASTTREE_EXEC, else PATH.

args <- commandArgs(trailingOnly = TRUE)
from_env <- function(var, name) {
  p <- Sys.getenv(var, "")
  if (nzchar(p)) p else unname(Sys.which(name))
}
raxml <- from_env("RAXML_NG_EXEC", "raxml-ng")
fasttree <- from_env("FASTTREE_EXEC", "FastTree")
run <- function(exec, a, ...) {
  status <- system2(exec, shQuote(a, type = "cmd"), ...)
  if (is.null(status)) 0L else as.integer(status)
}

modes <- c("--ancestral", "--evaluate", "--bootstrap", "--support", "--check", "--parse", "--rfdist")
if (any(args %in% modes)) quit(status = run(raxml, args))

# Dowser writes these as -msa and -prefix; accept both spellings.
value_of <- function(flags) {
  at <- which(args %in% flags)
  if (length(at) && at[1] < length(args)) args[at[1] + 1] else ""
}
msa <- value_of(c("-msa", "--msa"))
prefix <- value_of(c("-prefix", "--prefix"))
if (!nzchar(msa) || !nzchar(prefix)) {
  message("raxml-ng-fasttree: no --msa or --prefix in: ", paste(args, collapse = " "))
  quit(status = 2)
}

# FastTree's phylip reader is strict on name width and dowser writes relaxed phylip, so convert to FASTA.
draft_fa <- paste0(prefix, ".fasttree.fa")
draft <- paste0(prefix, ".fasttree.tree")
rows <- strsplit(trimws(readLines(msa)[-1]), "[[:space:]]+")
rows <- Filter(function(r) length(r) >= 2, rows)
writeLines(unlist(lapply(rows, function(r) c(paste0(">", r[1]), r[2]))), draft_fa)

# Same FastTree settings as the bash wrapper.
status <- run(fasttree, c("-nt", "-gtr", "-spr", "4", "-mlacc", "2", "-slownni",
                          "-nosupport", "-quiet", "-nopr", draft_fa), stdout = draft)
unlink(draft_fa)
if (status != 0L) quit(status = status)

quit(status = run(raxml, c(args, "--evaluate", "--tree", draft)))
