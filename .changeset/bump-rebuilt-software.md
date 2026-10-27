---
"@platforma-open/milaboratories.lineage-trees.immcantation": patch
"@platforma-open/milaboratories.lineage-trees.workflow": patch
"@platforma-open/milaboratories.lineage-trees.block": patch
---

Update MiXCR, MiTool to the versions with the rebuilt Docker images. The old images trust the compromised eMudhra intermediate CAs. Rebuild the immcantation Docker image, so it gets the current ca-certificates package.
