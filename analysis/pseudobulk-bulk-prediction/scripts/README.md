This directory contains scripts used in the analysis to predict bulk expression from matched pseudobulk expression.

* `sync-data-files.R`: This script syncs data files for analysis.
  * It identifies solid tumor samples of interest libraries of interest and obtains their bulk `quant.sf` files and single-cell `_processed.rds` files from S3.
  Files are saved in `../data/scpca-data/<project id>/<sample id>/`.
* `calculate-tpm.R`: This script creates a TSV of TPM values for all samples in a given project.
TSV files are exported to `../data/tpm/<project id>-tpm.tsv`.
* `calculate-pseudobulk.R`: This script creates a TSV of pseudobulk values, using two different approaches, for all samples in a given project.
TSV files are exported to `../data/pseudobulk/<project id>-pseudobulk.tsv`.
