---
'@platforma-open/milaboratories.lineage-trees.model': minor
'@platforma-open/milaboratories.lineage-trees.ui': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Collect tree nodes into named baskets.

Select rows in a tree's node table, or take a mutational path whole or in part,
and add them to a new or existing basket. Each basket is its own section: a node
table across lineages, with rows that can be removed, a name that can be edited,
and the table's export button.

Baskets are a view over this block's trees and reach no other block, so editing
one never stales the run. Each entry keeps a copy of the node's label and
reconstructed sequences. After a rerun on different settings, which numbers the
trees afresh, entries from the earlier run are listed from that copy rather than
matched to whatever node now carries their ids.
