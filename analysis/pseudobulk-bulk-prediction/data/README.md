This directory contains data files used as inputs to analysis.
Most contents of this directory are ignored by Git, but files included in Git are described here.

* `bulk-library-sample-ids.tsv`: a TSV of all bulk RNA-seq libraries and samples included in analysis.
  * This file was created with `../scripts/sync-data-files.R`
* `{project_id}_panglao-genesets.tsv` files contain gene sets for use in overrepresentation analysis.
  * Files were created by `../scripts/prepare-ora-gene-sets.R`
