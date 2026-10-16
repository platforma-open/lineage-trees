---
'@platforma-open/milaboratories.lineage-trees.block': patch
'@platforma-open/milaboratories.lineage-trees.model': patch
'@platforma-open/milaboratories.lineage-trees.ui': patch
'@platforma-open/milaboratories.lineage-trees.workflow': patch
'@platforma-open/milaboratories.lineage-trees.hilary': patch
'@platforma-open/milaboratories.lineage-trees.immcantation': patch
---

Rename the block to BCR Lineage Trees and mark the IgPhyML setting as slow. Donors without clonotypes in the picked datasets no longer run any stage, the donor progress bars move without a click, the sequence search counts its hits in a fraction of the memory, and the trees log drops two Dowser warnings that fired for every IgPhyML lineage.
