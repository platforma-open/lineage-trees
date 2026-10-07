#!/usr/bin/env Rscript
# Build MiXCR-shaped inputs for trees.R from dowser's paired example data.
# Each scenario packs many cases into one run, since starting R with dowser is the main cost.
# Scenarios under <out>/:
#   paired     AIRR route with light chains: split light V/J, light-less clones, known antibody, insertions, "/" in ids
#   joins      AIRR route with an ambiguous tuple join
#   table      heavy only, alignments on the clonotype table; deep lineage and exact copies
#   truncated  table route with partial 5' coverage
#   twins      one heavy chain with two light chains, both sharing one lineage
#   tiny       two tips per lineage, for IgPhyML
#   empty      a donor with no data
#   bare       no alignment columns and no AIRR export

suppressMessages(library(dplyr))
set.seed(7)

out_root <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_root)) stop("usage: fixture.R <output dir>")

# Load the data only; attaching dowser is slow.
data("ExampleMixedDb", package = "dowser", envir = environment())
db <- ExampleMixedDb
ungap <- function(s) toupper(gsub("[.-]", "", s))
key <- function(cell) paste0("ck_", cell)

# One clonotype per cell with distinct alignments, as MiXCR would give.
heavy <- db %>% filter(locus == "IGH") %>% distinct(cell_id, .keep_all = TRUE) %>%
  distinct(sequence_alignment, .keep_all = TRUE) %>% as.data.frame()
light <- db %>% filter(locus != "IGH") %>% distinct(cell_id, .keep_all = TRUE) %>% as.data.frame()

AIRR_COLS <- c("sequence_id", "cell_id", "v_call", "j_call", "junction",
               "sequence_alignment", "germline_alignment")

# `alignments`: "airr" writes AIRR exports, "table" puts them on the clonotype table, "none" neither.
write_scenario <- function(dir, heavy, light, light_less, with_light_columns = TRUE,
                           clone_prefix = "c", alignments = "airr") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  lh <- light[match(heavy$cell_id, light$cell_id), ]
  lh[heavy$cell_id %in% light_less, ] <- NA
  blank <- function(x) ifelse(is.na(x), "", x)

  clono <- data.frame(
    sequence_id = key(heavy$cell_id),
    v_call = heavy$v_call, j_call = heavy$j_call, junction = heavy$junction,
    main_sequence = ungap(heavy$sequence_alignment),
    stringsAsFactors = FALSE)
  if (alignments == "table") {
    clono$sequence_alignment <- heavy$sequence_alignment
    clono$germline_alignment <- heavy$germline_alignment
  }
  if ("is_known" %in% names(heavy)) clono$is_known <- heavy$is_known
  if (with_light_columns) {
    clono$v_call_light <- blank(lh$v_call)
    clono$j_call_light <- blank(lh$j_call)
    clono$junction_light <- blank(lh$junction)
    clono$main_sequence_light <- blank(ungap(lh$sequence_alignment))
  }
  write.table(clono, file.path(dir, "clonotypes.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(data.frame(sequence_id = key(heavy$cell_id),
                         clone_id = paste0(clone_prefix, heavy$clone_id)),
              file.path(dir, "clones.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (alignments == "airr") {
    dir.create(file.path(dir, "airr"), showWarnings = FALSE)
    # sequence_id is a per-file index and cell_id is a decoy; the join must use the gene/junction tuple.
    rows <- rbind(heavy[, AIRR_COLS], if (with_light_columns)
      light[!(light$cell_id %in% light_less) & light$cell_id %in% heavy$cell_id, AIRR_COLS]
      else light[0, AIRR_COLS])
    rows$sequence_id <- paste0("clone.", seq_len(nrow(rows)))
    half <- ceiling(nrow(rows) / 2)
    write.table(rows[1:half, ], file.path(dir, "airr", "s1.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(rows[(half + 1):nrow(rows), ], file.path(dir, "airr", "s2.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  writeLines(key(light_less), file.path(dir, "light-less.txt"))
}
note <- function(dir, name, ids) writeLines(ids, file.path(out_root, dir, name))

sizes <- sort(table(heavy$clone_id), decreasing = TRUE)
clone_of <- function(k) heavy$cell_id[heavy$clone_id == names(sizes)[k]]

# --- paired ---------------------------------------------------------------
split_clone <- clone_of(1)
recombined <- split_clone[seq_len(ceiling(length(split_clone) / 2))]
split_light <- light
split_light$v_call[split_light$cell_id %in% recombined] <- "IGLV9-49*01"
split_light$j_call[split_light$cell_id %in% recombined] <- "IGLJ3*02"
# Two clones with no light chain; the second drops a member to fall under RAxML's 4-sequence floor.
dark <- c(clone_of(2), clone_of(3)[-1])
known_clone <- clone_of(4)
insert_at <- function(aln, pos, what) {
  s <- strsplit(aln, "")[[1]]
  paste(c(s[1:pos], what, s[(pos + 1):length(s)]), collapse = "")
}
paired <- heavy[heavy$cell_id != clone_of(3)[1], ]
paired$is_known <- ifelse(paired$cell_id == known_clone[1], "true", "false")
# FR1 insertions ("-" in germline): 3nt off a codon boundary, 6nt on one. Plus one 1nt deletion in the query.
odd <- seq(1, nrow(paired), by = 2)
for (i in odd) {
  n <- if (i %% 4 == 1) 3 else 6
  pos <- if (n == 3) 31 else 30
  paired$sequence_alignment[i] <- insert_at(paired$sequence_alignment[i], pos, rep("A", n))
  paired$germline_alignment[i] <- insert_at(paired$germline_alignment[i], pos, rep("-", n))
}
substr(paired$sequence_alignment[odd[1]], 60, 60) <- "-"
write_scenario(file.path(out_root, "paired"), paired, split_light,
               c(dark, known_clone[2]), clone_prefix = "mouse D0_1/")
note("paired", "recombined.txt", key(recombined))
note("paired", "split-clone.txt", paste0("mouse D0_1/", names(sizes)[1]))
note("paired", "dark.txt", key(dark))
note("paired", "known.txt", key(known_clone[1]))
note("paired", "inserted-rows.txt", as.character(length(odd)))

# --- bulk -----------------------------------------------------------------
# The split clone plus a bulk copy of each member: no light chain, 120 bases short on
# the 5' side, one base changed so it joins its own AIRR row.
bulk_of <- heavy[heavy$cell_id %in% split_clone, ]
bulk_of$cell_id <- paste0("bulk_", bulk_of$cell_id)
for (col in c("sequence_alignment", "germline_alignment")) {
  bulk_of[[col]] <- substr(bulk_of[[col]], 121, nchar(bulk_of[[col]]))
}
flip <- c(A = "C", C = "A", G = "T", T = "G")
at <- vapply(strsplit(bulk_of$sequence_alignment, ""), function(s) which(s %in% names(flip))[5], integer(1))
substr(bulk_of$sequence_alignment, at, at) <- flip[substr(bulk_of$sequence_alignment, at, at)]
write_scenario(file.path(out_root, "bulk"), rbind(heavy[heavy$cell_id %in% split_clone, ], bulk_of),
               split_light, bulk_of$cell_id)
note("bulk", "recombined.txt", key(recombined))
note("bulk", "bulk.txt", key(bulk_of$cell_id))

# --- joins ----------------------------------------------------------------
# Clone members share one junction (as with assembly wider than CDR3), and one light chain pairs with all.
splice_junction <- function(aln, from, to) {
  if (nchar(from) != nchar(to)) return(aln)
  sub(from, to, aln, fixed = TRUE)
}
wide <- heavy %>% group_by(clone_id) %>%
  mutate(sequence_alignment = mapply(splice_junction, sequence_alignment, junction, junction[1]),
         junction = junction[1]) %>%
  ungroup() %>% as.data.frame()
one_light <- light[rep(1, nrow(heavy)), ]
one_light$cell_id <- heavy$cell_id
write_scenario(file.path(out_root, "joins"), wide, one_light, character())

# --- table ----------------------------------------------------------------
# Pad the biggest clone with mutants, and copy two members exactly.
# One copy is a known antibody; each pair must merge into one tip that keeps the known antibody and both links.
point_mutate <- function(aln, positions) {
  s <- strsplit(aln, "")[[1]]
  usable <- which(s %in% c("A", "C", "G", "T"))
  for (i in positions) {
    k <- usable[((i * 7) %% length(usable)) + 1]
    s[k] <- setdiff(c("A", "C", "G", "T"), s[k])[1]
  }
  paste(s, collapse = "")
}
deep_source <- heavy[heavy$clone_id == names(sizes)[1], ]
extra <- do.call(rbind, lapply(seq_len(8), function(i) {
  row <- deep_source[((i - 1) %% nrow(deep_source)) + 1, ]
  row$cell_id <- paste0(row$cell_id, "_deep", i)
  row$sequence_alignment <- point_mutate(row$sequence_alignment, seq_len(i))
  row
}))
deep <- rbind(heavy, extra)
dup <- deep_source[1:2, ]
dup$cell_id <- paste0(dup$cell_id, "_copy")
table_db <- rbind(deep, dup)
table_db$is_known <- ifelse(table_db$cell_id == dup$cell_id[1], "true", "false")
write_scenario(file.path(out_root, "table"), table_db, light, character(),
               with_light_columns = FALSE, alignments = "table")
note("table", "duplicated.txt", key(c(deep_source$cell_id[1:2], dup$cell_id)))
note("table", "known.txt", key(dup$cell_id[1]))

# --- twins ----------------------------------------------------------------
# The deep clone plus twins of three members: the same heavy chain, a light chain one FR1
# base apart. A twin is a distinct clonotype but not a distinct heavy sequence.
twin_of <- deep_source[deep_source$cell_id %in% light$cell_id, ][1:3, ]
twins <- twin_of
twins$cell_id <- paste0(twins$cell_id, "_twin")
twin_light <- light[match(twin_of$cell_id, light$cell_id), ]
twin_light$cell_id <- twins$cell_id
at <- vapply(strsplit(twin_light$sequence_alignment, ""), function(s) which(s %in% names(flip))[5], integer(1))
substr(twin_light$sequence_alignment, at, at) <- flip[substr(twin_light$sequence_alignment, at, at)]
write_scenario(file.path(out_root, "twins"), rbind(deep, twins), rbind(light, twin_light), character())
note("twins", "twins.txt", key(twins$cell_id))
note("twins", "sources.txt", key(twin_of$cell_id))

# --- truncated ------------------------------------------------------------
# Cuts of 115 and 170 leave the start off a codon boundary; a "." pad over the uncovered 5' side must be cut.
cut_left <- function(d, n) {
  for (col in c("sequence_alignment", "germline_alignment")) d[[col]] <- substr(d[[col]], n + 1, nchar(d[[col]]))
  d
}
truncated <- heavy
kind <- seq_len(nrow(truncated)) %% 4
truncated[kind == 1, ] <- cut_left(truncated[kind == 1, ], 115)
truncated[kind == 2, ] <- cut_left(truncated[kind == 2, ], 170)
pad <- kind == 3
truncated$sequence_alignment[pad] <- paste0(strrep(".", 115),
  substr(truncated$sequence_alignment[pad], 116, nchar(truncated$sequence_alignment[pad])))
write_scenario(file.path(out_root, "truncated"), truncated, light, character(),
               with_light_columns = FALSE, alignments = "table")

# --- marks -----------------------------------------------------------------
# The table scenario with known antibodies in a separate annotations table, as merge now writes them.
marks_dir <- file.path(out_root, "marks")
dir.create(marks_dir, showWarnings = FALSE)
file.copy(list.files(file.path(out_root, "table"), full.names = TRUE), marks_dir, recursive = TRUE)
marks_clono <- read.delim(file.path(marks_dir, "clonotypes.tsv"), colClasses = "character")
write.table(data.frame(sequence_id = marks_clono$sequence_id, data_source = "table",
                       is_known = marks_clono$is_known),
            file.path(marks_dir, "annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(marks_clono[, setdiff(names(marks_clono), "is_known")],
            file.path(marks_dir, "clonotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# --- oddids ----------------------------------------------------------------
# Ids holding the characters tree files rewrite: ":", ";", ",", "=" and a space.
odd <- function(d) { d$cell_id <- paste0("k:a;b,c=d e ", d$cell_id); d }
write_scenario(file.path(out_root, "oddids"), odd(heavy), odd(light), character())

# --- gapped, mixedgaps -----------------------------------------------------
# One donor's rows IMGT-gapped, then half of them ungapped as another tool would write them.
write_scenario(file.path(out_root, "gapped"), heavy, light, character(),
               with_light_columns = FALSE, alignments = "table")
ungap_row <- function(s, g) {
  sc <- strsplit(s, "")[[1]]; gc <- strsplit(g, "")[[1]]
  keep <- gc != "."
  c(paste(sc[keep], collapse = ""), paste(gc[keep], collapse = ""))
}
mixedgaps <- heavy
half <- seq_len(nrow(mixedgaps)) %% 2 == 0
pairs <- mapply(ungap_row, mixedgaps$sequence_alignment[half], mixedgaps$germline_alignment[half],
                USE.NAMES = FALSE)
mixedgaps$sequence_alignment[half] <- pairs[1, ]
mixedgaps$germline_alignment[half] <- pairs[2, ]
write_scenario(file.path(out_root, "mixedgaps"), mixedgaps, light, character(),
               with_light_columns = FALSE, alignments = "table")

# --- tiny, empty, bare ----------------------------------------------------
tiny <- heavy %>% group_by(clone_id) %>% slice_head(n = 2) %>% ungroup() %>% as.data.frame()
write_scenario(file.path(out_root, "tiny"), tiny, light, character())

empty_dir <- file.path(out_root, "empty")
dir.create(file.path(empty_dir, "airr"), recursive = TRUE, showWarnings = FALSE)
writeLines("sequence_id\tv_call\tj_call\tjunction\tmain_sequence", file.path(empty_dir, "clonotypes.tsv"))
writeLines("sequence_id\tclone_id", file.path(empty_dir, "clones.tsv"))
writeLines(paste(AIRR_COLS, collapse = "\t"), file.path(empty_dir, "airr", "s1.tsv"))

write_scenario(file.path(out_root, "bare"), heavy, light, character(),
               with_light_columns = FALSE, alignments = "none")
