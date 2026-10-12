---
'@platforma-open/milaboratories.lineage-trees.workflow': patch
'@platforma-open/milaboratories.lineage-trees.block': patch
---

Fix exportAirr running out of memory on large samples. Memory now grows with the clns size, from 8 GiB.
