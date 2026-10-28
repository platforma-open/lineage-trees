# @platforma-open/milaboratories.lineage-trees.block

## 1.3.0

### Minor Changes

- 20f6475: Accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) per input dataset. Each picked dataset that has subset columns gets its own optional "Subset of …" picker, so one dataset can come in whole and another narrowed to a subset.

  A subset restricts its dataset strictly, before anything is computed: only its clonotypes enter clustering, allele inference, alignment and tree building, so lineages are never linked through a clonotype outside the subset. MiXCR alignments are still read from the whole clns, but only rows matching an included clonotype are used.

  Unlike the single-input blocks, the subset is chosen per dataset. Each dataset's exported columns (clonotype-to-lineage links, germline mutations, known antibody distances and links, and baskets) carry `pl7.app/inputSubset` and a trace rooted in that dataset's subset; a known antibody link names the candidate dataset's subset. Baskets can hold nodes saved in earlier runs, so on a subset run each exported basket is cut to its dataset's subset, and the basket page hides nodes from earlier runs while a subset is selected. Lineage-level columns are not stamped, since a lineage can combine several datasets with different subsets. Abundance totals and fractions of a subset run are relative to the clonotypes that entered, so they differ from a full-data run. Runs without a subset are unchanged.

## 1.2.4

### Patch Changes

- c2db88d: Update MiXCR, MiTool to the versions with the rebuilt Docker images. The old images trust the compromised eMudhra intermediate CAs. Rebuild the immcantation Docker image, so it gets the current ca-certificates package.

## 1.2.3

### Patch Changes

- a9cdf0e: Name the known antibody link "Nearest known antibody", so downstream blocks show "Reached via Nearest known antibody" on the link icon. Update to SDK 1.86.0.

## 1.2.2

### Patch Changes

- cd9c269: update sdk dependencies

## 1.2.1

### Patch Changes

- a6e8755: Fix exportAirr running out of memory on large samples. Memory now grows with the clns size, from 8 GiB.

## 1.2.0

### Minor Changes

- 863caf4: Add isotypes, mutations from the germline and the MRCA, and reads, UMIs and cells to the trees.

  Each observed node carries an isotype, from the producer's isotype column or, for imported AIRR
  data, the class of its heavy C gene. Where the node's clonotypes differ in isotype, the one with
  the most cells, UMIs or reads wins, in the first of those every clonotype has; with none shared,
  the node is left blank rather than comparing counts in different units. The tree tooltip shows
  it when a run has isotypes.

  Every node also gets its changes since the germline and since the most recent common ancestor
  of the lineage's observed sequences, per chain, in nucleotides and amino acids: the list, the
  count and, from the germline, the rate over positions settled at both ends. Changes from the
  germline are in V and J only, as the germline has no junction; MRCA figures are blank above it.
  Nodes also get their distance from the germline along the tree and their heavy chain V and J
  identity to the germline.

  Reads, UMIs and cells are summed per node and per lineage, each from the datasets that count
  it, and blank where none does. The node's primary abundance is now shown only when every
  dataset counts in the same unit, rather than adding reads to UMIs.

  The tree's node table now lists observed sequences only and opens on their V, D and J alleles,
  CDR3 and VDJRegion. Steps and mutation columns, including those from the MRCA, open on the
  mutational path and basket tables instead.

## 1.1.0

### Minor Changes

- fc76d58: Remove the Heavy AA consensus distance (aa-cdist) column and the Consensus sequences
  count that explained it.

  The metric reproduces its authors' published performance on their own simulated
  families, but we could not reproduce it on real data: it is anti-predictive on the
  Landais PCT64 lineage using the authors' deposited sequences, null across 119
  germinal centres in the gcreplay dataset, and the second real-data example in the
  source paper cannot be reproduced at all because the clonal family composition was
  never published. The metric assumes a family converging on a single optimum in which
  mutations are net deleterious; neither condition holds in the real lineages tested,
  and the block cannot detect when they fail. Clonotypes keep their germline mutation
  count and, with anchors, their mutations to the nearest anchor.
