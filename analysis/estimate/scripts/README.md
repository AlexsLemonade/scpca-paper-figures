This directory contains scripts used in the analysis to predict bulk expression from matched pseudobulk expression.

* `sync-data-files.R`: This script syncs data files for analysis.
  * It identifies solid tumor samples of interest libraries of interest and obtains their bulk `quant.sf` files and single-cell `_processed.rds` files from S3.
  Files are saved in `../data/scpca/<project id>/<sample id>/`.
  * It also exports two TSVs to `../data`:
    * Map for bulk library and sample ids
    * Map for ensembl ids and gene symbols
