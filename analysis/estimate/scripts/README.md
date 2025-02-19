This directory contains scripts used in the analysis to predict bulk expression from matched pseudobulk expression.

* `sync-data-files.R`: This script syncs data files for analysis.
  * It identifies solid tumor libraries of interest and obtains their bulk `quant.sf` files and single-cell `_processed.rds` files from S3.
  Files are saved in `../data/scpca/<project id>/<sample id>/`.
  * It also exports a TSV to `../data` that maps ensembl ids to gene symbols
* `prepare-expression-data.R`: This script prepares data frames of expression for a single project for input to `ESTIMATE`
  * The script exports an RDS file with a list of data frames for bulk counts, bulk TPM, and pseudobulk counts to `../data`
