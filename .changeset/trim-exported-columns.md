---
'@platforma-open/milaboratories.lineage-trees.workflow': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Export fewer anchor and ranking columns.

Anchor sets no longer get Heavy AA consensus distance or any anchor column, all
blank or meaningless for the reference itself. Bulk datasets no longer get light
chain anchor columns, which were always blank. The NT mutations to anchor are no
longer emitted: a silent change says nothing about binding.

The nearest-anchor linker now carries `pl7.app/linkLabel: "Anchor"`, so a
consumer that names columns by the linker that reached them calls the anchor's
columns "Anchor" rather than "Nearest anchor".
