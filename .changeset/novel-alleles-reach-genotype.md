---
'@platforma-open/milaboratories.lineage-trees.immcantation': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Use the novel alleles TIgGER finds. Until now none of them reached the genotype.

TIgGER's `inferGenotype` adds novel alleles to the candidate set only when
`find_unmutated = TRUE`. With the flag off, the genotype held only alleles
already in `v_call`, which on a MiXCR `*00` library is one per gene. On 16
rhesus macaque datasets the stage found 50 novel alleles, used none of them,
and rewrote 0 rows.

The flag is back on. The pruning it brings is handled by the guards added
since it was turned off: a row whose gene was pruned keeps its call, and a row
that would end further from its germline keeps its call. On one rhesus IgM
sample this takes total V mutations to 65% of baseline, with no row increased.
