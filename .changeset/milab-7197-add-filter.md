---
"@platforma-open/milaboratories.lineage-trees.block": minor
"@platforma-open/milaboratories.lineage-trees.model": minor
"@platforma-open/milaboratories.lineage-trees.ui": minor
"@platforma-open/milaboratories.lineage-trees.workflow": minor
---

Accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) per input dataset. Each picked dataset that has subset columns gets its own optional "Subset of …" picker, so one dataset can come in whole and another narrowed to a subset.

A subset restricts its dataset strictly, before anything is computed: only its clonotypes enter clustering, allele inference, alignment and tree building, so lineages are never linked through a clonotype outside the subset. MiXCR alignments are still read from the whole clns, but only rows matching an included clonotype are used.

Unlike the single-input blocks, the subset is chosen per dataset. Each dataset's exported columns (clonotype-to-lineage and node links, germline mutations, known antibody distances and links, member counts) carry `pl7.app/inputSubset` and a trace rooted in that dataset's subset. Lineage-level columns are not stamped, since a lineage can combine several datasets with different subsets. Runs without a subset are unchanged.
