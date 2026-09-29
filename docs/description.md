# Overview

Lineage Trees groups B cell clonotypes into clonal families and builds a family tree for each one. Members of a family come from the same starting B cell and differ by mutations picked up during an immune response. The block shows how each family developed and gives each clonotype numbers you can rank candidates by.

It works on antibody (IG) data only, bulk heavy chain or paired single cell, and can mix both in one run. Input comes from MiXCR Clonotyping or Import V(D)J Data. Families never cross donors: pick a sample metadata column that names each sample's donor, or all samples are treated as one donor.

What it does:

- **Families.** Clonotypes are grouped by their heavy chain with HILARy. Where light chains are present, families are then split by light chain.
- **Germline alleles.** TIgGER looks for V gene alleles the reference does not have. A sequence moves to a new allele only when it fits that allele better.
- **Trees.** Each family gets a tree with its ancestral sequences reconstructed, using FastTree and RAxML-NG through Dowser. IgPhyML is optional.
- **Anchors.** Mark a dataset of known antibodies as anchors, and every other member of an anchor's family gets its number of amino acid changes to the nearest anchor, counted along the tree.

The main ranking column, **Heavy AA consensus distance**, counts how many amino acids in a clonotype's heavy chain differ from its family's consensus. Lower is better. It is only given for families with at least 10 distinct heavy sequences, because a smaller family has no reliable consensus.

The ranking columns and family ids go to downstream blocks such as Lead Selection, which can rank candidates by them.

The consensus distance follows Ralph and Matsen (2020), *PLOS Computational Biology* 16(11): e1008391. Please also cite HILARy, TIgGER, Dowser, FastTree, RAxML-NG and, if used, IgPhyML; the block's README lists the references.
