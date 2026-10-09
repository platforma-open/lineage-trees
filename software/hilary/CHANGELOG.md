# @platforma-open/milaboratories.lineage-trees.hilary

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
