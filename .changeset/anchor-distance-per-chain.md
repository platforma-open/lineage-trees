---
'@platforma-open/milaboratories.lineage-trees.workflow': minor
'@platforma-open/milaboratories.lineage-trees.model': minor
'@platforma-open/milaboratories.lineage-trees.ui': minor
'@platforma-open/milaboratories.lineage-trees.immcantation': minor
'@platforma-open/milaboratories.lineage-trees.hilary': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Count mutations per chain and in amino acids, and rank on amino acids.

A paired tree is built over heavy and light joined end to end, and every
mutation count was taken over the join. A lineage with light chains reported
heavy plus light, one without reported heavy alone, under one label; and a
member without a light chain looked closer to an anchor than a paired sibling
just as far away, because its missing light half settles nothing.

Every branch is now counted per chain, split where dowser joined them, and each
chain is translated on its own so its frame does not rest on the other. The
node table carries heavy and light amino acid mutations, their counts, and the
nucleotide mutations, counts and unresolved positions beside them. Light
columns exist only when light chains were used.

The anchor distance becomes one column per chain, both to the one nearest
anchor, now chosen by the heavy chain's amino acid changes: Heavy AA mutations
to anchor and, on paired datasets, Light AA mutations to anchor are ranking
columns. A light figure is blank unless both ends carry a light chain. Heavy is always named and the chain is in the
domain, so a column reads the same in every run and never passes for a
two-chain total.

The reconstructed sequence is split the same way, into Heavy reconstructed
sequence and Light reconstructed sequence, the light one only when light chains
were used and blank in lineages without them. Both close every table, heavy
first, and the tree tooltip shows both.

aa-cdist, heavy chain only as before, is relabelled Heavy AA consensus distance.
