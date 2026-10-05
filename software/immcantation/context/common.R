# Helpers shared by trees.R and alleles.R, sourced from beside them.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, required = TRUE) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    if (required) stop("missing ", flag)
    return(NULL)
  }
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

read_tsv <- function(p) read.delim(p, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
# A call's gene: its first hit, allele stripped.
gene_of <- function(x) sub("\\*.*$", "", sub(",.*$", "", x))

# Elapsed wall time as H:MM:SS, from the start of this run unless told otherwise.
STARTED <- Sys.time()
clock <- function(since = STARTED) {
  s <- as.integer(round(as.numeric(difftime(Sys.time(), since, units = "secs"))))
  sprintf("%d:%02d:%02d", s %/% 3600L, s %/% 60L %% 60L, s %% 60L)
}

# Progress prefix; the block shows the last such line of a step's output.
PROGRESS_PREFIX <- "[==PROGRESS==]"
# Stamped with the elapsed time up front, since the UI reads the percentage off the end.
# Flushed, so the block sees the step as it starts.
progress <- function(text) { cat(sprintf("%s [%s] %s\n", PROGRESS_PREFIX, clock(), text)); flush(stdout()) }
