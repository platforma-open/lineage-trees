---
'@platforma-open/milaboratories.lineage-trees.workflow': minor
'@platforma-open/milaboratories.lineage-trees.model': minor
'@platforma-open/milaboratories.lineage-trees.ui': minor
'@platforma-open/milaboratories.lineage-trees.hilary': minor
'@platforma-open/milaboratories.lineage-trees.immcantation': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Add isotypes, mutations from the germline and the MRCA, and reads, UMIs and cells to the trees.

Each observed node carries an isotype, from the producer's isotype column or, for imported AIRR
data, the class of its heavy C gene. Where the node's clonotypes differ in isotype, the one with
the most cells, UMIs or reads wins, in the first of those every clonotype has; with none shared,
the node is left blank rather than comparing counts in different units. The tree tooltip shows
it when a run has isotypes.

Every node also gets its changes since the germline and since the most recent common ancestor
of the lineage's observed sequences, per chain, in nucleotides and amino acids: the list, the
count and, from the germline, the rate over positions settled at both ends. Changes from the
germline are in V and J only, as the germline has no junction; MRCA figures are blank above it.
Nodes also get their distance from the germline along the tree and their heavy chain V and J
identity to the germline.

Reads, UMIs and cells are summed per node and per lineage, each from the datasets that count
it, and blank where none does. The node's primary abundance is now shown only when every
dataset counts in the same unit, rather than adding reads to UMIs.

The tree's node table now lists observed sequences only and opens on their V, D and J alleles,
CDR3 and VDJRegion. Steps and mutation columns, including those from the MRCA, open on the
mutational path and basket tables instead.
