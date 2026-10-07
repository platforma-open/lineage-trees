# Overview

Groups B cell clonotypes into clonal lineages and builds a tree for each lineage, showing how it developed during an immune response and where each clonotype sits in that history. Members of a lineage descend from the same B cell and differ by the mutations they picked up along the way. The block works on antibody (IG) data, bulk heavy chain or paired single cell, from MiXCR Clonotyping or Import V(D)J Data, and can mix both in one run. Lineages never cross donors: a sample metadata column can name each sample's donor, otherwise all samples are treated as one donor.

Clonotypes are grouped into lineages by their heavy chain with HILARy, using a fixed similarity threshold or its adaptive mode, and split by light chain where light chains are present. TIgGER first infers the donor's V gene alleles. Each lineage then gets a tree with its ancestral sequences reconstructed, built with FastTree and RAxML-NG through Dowser, or optionally with IgPhyML.

Every clonotype gets its heavy chain mutation count from the germline. Datasets of known antibodies can be marked as such, and every other member of their lineage gets its distance to the nearest known antibody. The ranking columns and lineage ids go to downstream blocks such as Lead Selection.

When using this block in your research, cite the publications for the tools your run relied on, listed below.

> Spisak, N., Athènes, G., Dupic, T., Mora, T., & Walczak, A. M. (2024). Combining mutation and recombination statistics to infer clonal families in antibody repertoires. _eLife_ **13**, e86181. [https://doi.org/10.7554/eLife.86181](https://doi.org/10.7554/eLife.86181)

> Gadala-Maria, D., Yaari, G., Uduman, M., & Kleinstein, S. H. (2015). Automated analysis of high-throughput B-cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles. _PNAS_ **112**(8), E862-E870. [https://doi.org/10.1073/pnas.1417683112](https://doi.org/10.1073/pnas.1417683112)

> Hoehn, K. B., Pybus, O. G., & Kleinstein, S. H. (2022). Phylogenetic analysis of migration, differentiation, and class switching in B cells. _PLOS Computational Biology_ **18**(4), e1009885. [https://doi.org/10.1371/journal.pcbi.1009885](https://doi.org/10.1371/journal.pcbi.1009885)

> Price, M. N., Dehal, P. S., & Arkin, A. P. (2010). FastTree 2: approximately maximum-likelihood trees for large alignments. _PLOS ONE_ **5**(3), e9490. [https://doi.org/10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)

> Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. _Bioinformatics_ **35**(21), 4453-4455. [https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

> Hoehn, K. B., Vander Heiden, J. A., Zhou, J. Q., Lunter, G., Pybus, O. G., & Kleinstein, S. H. (2019). Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination. _PNAS_ **116**(45), 22664-22672. [https://doi.org/10.1073/pnas.1906020116](https://doi.org/10.1073/pnas.1906020116)
