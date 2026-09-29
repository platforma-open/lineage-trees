---
'@platforma-open/milaboratories.lineage-trees.immcantation': patch
'@platforma-open/milaboratories.lineage-trees.block': patch
---

Build trees on Windows. The raxml-ng wrapper that serves the topology search from FastTree gets an R port, run there through a batch file since Windows has no bash, and dowser's own forked parallelism runs as one process on Windows only. Linux and macOS are unchanged.
