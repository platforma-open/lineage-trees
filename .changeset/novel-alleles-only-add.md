---
'@platforma-open/milaboratories.lineage-trees.immcantation': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Novel alleles only add to the reference; they never replace it.

- Every reference allele the genotype dropped is put back, so no gene is lost
  and a sequence keeps its call unless a same-gene allele sits strictly closer.
- Moves to another gene are refused.
- The novel allele scan stops 10 positions before the end of the V, and uses
  the data's own germline width rather than TIgGER's fixed 312.
