---
'@platforma-open/milaboratories.lineage-trees.immcantation': patch
'@platforma-open/milaboratories.lineage-trees.block': patch
---

Report what the genotype did with the novel alleles.

The allele log now lists, per gene, the alleles the genotype kept with their
sequence counts, the novel alleles it left out, and the genes it dropped. The
genotype and final summary lines say how many of the novel alleles found made
it into the genotype, and how many sequences took a new allele.

The closer-germline check now judges only rows offered a different allele.
Rows given back the allele they already had were counted as refused
reassignments, which made a run that applied nothing report identical mutation
totals on both sides.
