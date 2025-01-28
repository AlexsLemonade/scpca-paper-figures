This directory contains scripts used in the analysis to predict bulk expression from matched pseudobulk expression.

* `prepare-data-files.R`: This script prepares some input data for analysis:
  * It identifies solid tumor bulk libraries of interest and obtains their `quant.sf` files from S3.
  These files will be used to calculate TPM values as input to prediction.
    * Files are saved in `../data/salmon-quant-files/<project id>/<sample id>/<library id>`.
* `calculate-tpm.R`: This script creates a TSV of TPM values for all samples in a given project.
TSV files are exported to `../data/tpm/<project id>-tpm.tsv`.
