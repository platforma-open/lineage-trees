# Lineage Trees

Reconstruct BCR clonal lineages and find the relatives of a known antibody. This Platforma block groups clonotypes into clonal lineages per donor, builds a phylogenetic tree for each lineage with its ancestral sequences reconstructed, and reports where each clonotype sits in that history: its mutations from the germline, and how far it lies along the tree from a known antibody.

Open-source analysis block for Platforma, the biologics discovery platform by MiLaboratories. For the full no-code workflow, see [platforma.bio](https://platforma.bio/).

## What it does

Clonotypes descended from one naive B cell form a clonal lineage whose members differ by somatic hypermutation, and a candidate's worth often depends on its place in that lineage rather than on its sequence alone. This block recovers the lineages, their trees, and the numbers those trees make available.

**Lineages.** Clonotypes are clustered with HILARy inside each donor, since shared ancestry only exists within one immune system. Clustering runs on heavy chains, which lets one lineage hold bulk and single-cell members at once. Light chains are applied afterwards: wherever any picked dataset carries them, known antibody datasets included, Dowser's `resolveLightChains` splits each heavy-chain clone by light V, J and junction length, and places every heavy-only member with its nearest paired member.

**Germline alleles.** An allele the reference does not know reads as a mutation shared by a whole lineage. Each donor's V alleles are inferred with TIgGER from that donor's own sequences, pooled across its datasets, and the heavy chains' V calls and germlines are updated to them before clustering and tree building. A donor with too few sequences keeps the reference alleles.

**Trees.** Each lineage's tree is built inside Dowser. FastTree drafts the topology, and raxml-ng fits the substitution model, the branch lengths and the ancestral sequences onto it. A full topology search does not scale to real lineage sizes, so only the draft is FastTree's. IgPhyML's HLP19 codon model, which knows the hotspot biases of somatic hypermutation, is optional. Lineages under three tips, and builds that fail, fall back to maximum parsimony; the Tree builder column says which ran.

**Known antibodies.** Any picked dataset can be marked as a known antibody dataset: antibodies whose relatives you are looking for, such as characterised leads or hits from a screen. Every other member of a lineage holding one gets the amino acid mutations along the tree to the nearest known antibody, per chain, and the id of that known antibody.

## Inputs & outputs

* **Input:** one or more IG clonotyping datasets, bulk heavy chain or paired single cell, from [MiXCR Clonotyping](https://github.com/platforma-open/mixcr-clonotyping) or [Import V(D)J Data](https://github.com/platforma-open/import-vdj-data), clustered together. Optionally a sample metadata column naming each sample's donor.
* **Output:** per clonotype, its lineage, its heavy chain nucleotide mutations from the germline and, with known antibodies, its mutations to the nearest known antibody. Per lineage, size, genes, representative CDR3, abundance, known antibodies held, the tree builder and, when several datasets were picked, members per dataset. Per node, the reconstructed heavy and (on paired data) light sequences, the mutations acquired on the branch reaching it, and its place in the tree: depth, terminal branch fraction and the size of its parent's clade. Lineages, germline mutation counts, ranking columns and baskets are available to downstream blocks.

## Specifications

| | |
|---|---|
| Block title in app | Lineage Trees |
| Receptor | BCR (IG) only |
| Modalities | Bulk heavy chain and paired single cell, mixed in one run |
| Tools | [HILARy](https://github.com/statbiophys/HILARy), [TIgGER](https://tigger.readthedocs.io/), [Dowser](https://dowser.readthedocs.io/), [FastTree](https://morgannprice.github.io/fasttree/), [raxml-ng](https://github.com/amkozlov/raxml-ng), optionally [IgPhyML](https://github.com/immcantation/igphyml) |
| Alignments | MiXCR's own where a dataset carries a `clns`, otherwise the dataset's `sequence_alignment` and `germline_alignment`. Nothing is realigned |
| Ranking columns | Heavy AA mutations to known antibody, Light AA mutations to known antibody on paired data (lower is better) |
| Views | Donor overview with per-stage progress and logs, lineage table, dendrogram and node table per lineage, mutational path, baskets, sequence search |

## Settings

Everything below Donor column lives under **Advanced**.

| Setting | Default | Effect |
|---|---|---|
| Datasets | none | One or more IG datasets, all clustered together |
| Known antibody datasets | none | Which picked datasets hold antibodies whose relatives you want. Needs two or more datasets |
| Donor column | Single donor | Sample metadata naming each sample's subject. No lineage spans donors |
| Clustering | Fixed threshold | Fixed threshold, or HILARy's adaptive mode |
| Clustering threshold | `0.2` | Fraction of CDR3 length below which two clonotypes are single-linked. Fixed mode |
| Precision | `0.99` | Share of pairs placed together that HILARy aims to have truly related. Adaptive mode |
| Sensitivity | `0.9` | Share of truly related pairs that HILARy aims to place together. Adaptive mode |
| IgPhyML | None | Which lineages IgPhyML builds: none, those holding a known antibody (offered when known antibodies are picked), or all |
| Maximum tips per tree | unset | Subsamples larger lineages at random before building, always keeping known antibodies |
| Minimum tips per tree | unset | Smaller lineages keep their membership and get no tree |

A clonotype left out of a tree keeps its lineage and its germline mutation count, and loses only the tree columns.

## Use cases

* **Find relatives of a known antibody:** mark its dataset as a known antibody dataset and rank by Heavy AA mutations to known antibody.
* **Read a maturation history:** follow the mutational path from the germline to a candidate, and collect nodes of interest into baskets.
* **Join bulk and single-cell data:** cluster both for one donor into shared lineages, with light chains from the paired data.
* **Feed lead selection:** supply lineages and ranking columns to [Lead Selection](https://github.com/platforma-open/antibody-tcr-lead-selection).

## FAQ

### Fixed threshold or adaptive?

The fixed threshold single-links clonotypes on CDR3 within a V, J and CDR3-length class: fast and predictable. Adaptive infers a threshold per class for the precision and sensitivity you ask for, and adds a test on mutations shared outside the CDR3. It is calibrated on human repertoires, and the block warns when it is picked. In adaptive mode, a donor where any clonotype lacks an alignment is clustered on CDR3 alone, without the shared-mutation test.

### Why mutations along the tree, and not sequence identity?

Because the question is shared history. The distance along the reconstructed tree, through the most recent common ancestor, counts the changes that actually separate two clonotypes. Positions the reconstruction could not settle are left out rather than guessed, and silent changes are not counted, since they say nothing about binding.

### Are imported datasets realigned?

No. Realigning would give a dataset a second annotation that other tables in the project cannot match. Imported and MiXCR datasets annotated differently may therefore not cluster together, and after a run the block warns when common genes of one source never appear in the other.

## Citation

If you use this block in your research, please cite the tools your run relied on:

> **HILARy.** Spisak, N., Athènes, G., Dupic, T., Mora, T., & Walczak, A. M. (2024). Combining mutation and recombination statistics to infer clonal families in antibody repertoires. *eLife* **13**, e86181. [doi:10.7554/eLife.86181](https://doi.org/10.7554/eLife.86181)

> **Dowser.** Hoehn, K. B., Pybus, O. G., & Kleinstein, S. H. (2022). Phylogenetic analysis of migration, differentiation, and class switching in B cells. *PLOS Computational Biology* **18**(4), e1009885. [doi:10.1371/journal.pcbi.1009885](https://doi.org/10.1371/journal.pcbi.1009885)

> **FastTree 2.** Price, M. N., Dehal, P. S., & Arkin, A. P. (2010). FastTree 2: approximately maximum-likelihood trees for large alignments. *PLOS ONE* **5**(3), e9490. [doi:10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)

> **RAxML-NG.** Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics* **35**(21), 4453-4455. [doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

> **IgPhyML.** Hoehn, K. B., Vander Heiden, J. A., Zhou, J. Q., Lunter, G., Pybus, O. G., & Kleinstein, S. H. (2019). Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination. *PNAS* **116**(45), 22664-22672. [doi:10.1073/pnas.1906020116](https://doi.org/10.1073/pnas.1906020116)

> **TIgGER.** Gadala-Maria, D., Yaari, G., Uduman, M., & Kleinstein, S. H. (2015). Automated analysis of high-throughput B-cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles. *PNAS* **112**(8), E862-E870. [doi:10.1073/pnas.1417683112](https://doi.org/10.1073/pnas.1417683112)

## Part of the Platforma ecosystem

This block is part of [Platforma](https://platforma.bio/) by [MiLaboratories](https://github.com/milaboratory). Explore the other open-source blocks at [github.com/platforma-open](https://github.com/platforma-open) and the antibody discovery docs at [docs.platforma.bio](https://docs.platforma.bio/biology-guides/antibody-discovery/).
