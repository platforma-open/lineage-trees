---
'@platforma-open/milaboratories.lineage-trees.immcantation': minor
'@platforma-open/milaboratories.lineage-trees.block': minor
---

Run germline reconstruction, allele inference and tree building without Docker.

Both immcantation entrypoints now also ship as an R package on the
`runenv-r-lineage-trees` run environment, which also provides raxml-ng, FastTree
and IgPhyML, so the block runs on a backend that has no Docker. Where
Docker is available the existing image is still used.
