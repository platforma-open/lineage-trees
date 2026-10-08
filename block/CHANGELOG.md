# @platforma-open/milaboratories.lineage-trees.block

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
