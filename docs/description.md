# Overview

Groups B cell clonotypes into clonal families and builds a family tree for each one, showing how each family developed during an immune response and giving every clonotype numbers that candidates can be ranked by. Members of a family descend from the same B cell and differ by the mutations they picked up along the way. The block works on antibody (IG) data, bulk heavy chain or paired single cell, from MiXCR Clonotyping or Import V(D)J Data, and can mix both in one run. Families never cross donors: a sample metadata column can name each sample's donor, otherwise all samples are treated as one donor.

Clonotypes are grouped into families by their heavy chain with HILARy v1.2.4, using either a fixed similarity threshold or HILARy's adaptive mode, which sets its own thresholds from the data. Where light chains are present, families are then split by light chain. Before clustering, TIgGER v1.1.3 looks for V gene alleles that the reference lacks, and a sequence moves to a new allele only when it fits that allele better. Each family then gets a tree with its ancestral sequences reconstructed, built with FastTree v2.1.11 and RAxML-NG v2.0.3 through Dowser v2.5.1, or optionally with IgPhyML.

The main ranking column, Heavy AA consensus distance, counts the amino acids in a clonotype's heavy chain that differ from its family's consensus, following Ralph and Matsen (2020); lower is better. It is given only for families with at least 10 distinct heavy sequences, since a smaller family has no reliable consensus. A dataset of known antibodies can be marked as anchors, and every other member of an anchor's family then gets its number of amino acid changes to the nearest anchor, counted along the tree. The ranking columns and family ids go to downstream blocks such as Lead Selection.

When using this block in your research, cite the publications for the tools your run relied on, listed below.

> Spisak, N., Athènes, G., Dupic, T., Mora, T., & Walczak, A. M. (2024). Combining mutation and recombination statistics to infer clonal families in antibody repertoires. _eLife_ **13**, e86181. [https://doi.org/10.7554/eLife.86181](https://doi.org/10.7554/eLife.86181)

> Gadala-Maria, D., Yaari, G., Uduman, M., & Kleinstein, S. H. (2015). Automated analysis of high-throughput B-cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles. _PNAS_ **112**(8), E862-E870. [https://doi.org/10.1073/pnas.1417683112](https://doi.org/10.1073/pnas.1417683112)

> Hoehn, K. B., Pybus, O. G., & Kleinstein, S. H. (2022). Phylogenetic analysis of migration, differentiation, and class switching in B cells. _PLOS Computational Biology_ **18**(4), e1009885. [https://doi.org/10.1371/journal.pcbi.1009885](https://doi.org/10.1371/journal.pcbi.1009885)

> Price, M. N., Dehal, P. S., & Arkin, A. P. (2010). FastTree 2: approximately maximum-likelihood trees for large alignments. _PLOS ONE_ **5**(3), e9490. [https://doi.org/10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)

> Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. _Bioinformatics_ **35**(21), 4453-4455. [https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

> Hoehn, K. B., Vander Heiden, J. A., Zhou, J. Q., Lunter, G., Pybus, O. G., & Kleinstein, S. H. (2019). Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination. _PNAS_ **116**(45), 22664-22672. [https://doi.org/10.1073/pnas.1906020116](https://doi.org/10.1073/pnas.1906020116)

> Ralph, D. K., & Matsen IV, F. A. (2020). Using B cell receptor lineage structures to predict affinity. _PLOS Computational Biology_ **16**(11), e1008391. [https://doi.org/10.1371/journal.pcbi.1008391](https://doi.org/10.1371/journal.pcbi.1008391)
