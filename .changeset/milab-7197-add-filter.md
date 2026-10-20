---
"@platforma-open/milaboratories.lineage-trees.block": minor
"@platforma-open/milaboratories.lineage-trees.model": minor
"@platforma-open/milaboratories.lineage-trees.ui": minor
"@platforma-open/milaboratories.lineage-trees.workflow": minor
---

Accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) per input dataset. Each picked dataset that has subset columns gets its own optional "Subset of …" picker, so one dataset can come in whole and another narrowed to a subset.

A subset restricts its dataset strictly, before anything is computed: only its clonotypes enter clustering, allele inference, alignment and tree building, so lineages are never linked through a clonotype outside the subset. MiXCR alignments are still read from the whole clns, but only rows matching an included clonotype are used.

Unlike the single-input blocks, the subset is chosen per dataset. Each dataset's exported columns (clonotype-to-lineage links, germline mutations, known antibody distances and links, and baskets) carry `pl7.app/inputSubset` and a trace rooted in that dataset's subset; a known antibody link names the candidate dataset's subset. Baskets can hold nodes saved in earlier runs, so on a subset run each exported basket is cut to its dataset's subset, and the basket page hides nodes from earlier runs while a subset is selected. Lineage-level columns are not stamped, since a lineage can combine several datasets with different subsets. Abundance totals and fractions of a subset run are relative to the clonotypes that entered, so they differ from a full-data run. Runs without a subset are unchanged.
