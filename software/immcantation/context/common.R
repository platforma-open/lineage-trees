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
