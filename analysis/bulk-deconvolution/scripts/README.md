This directory contains scripts used in the bulk deconvolution analysis.

- `prepare-data-files.R`: This script prepares some input data for analysis:
  - First, it uses an ScPCA SCE to prepare and export a table mapping ensembl ids and gene symbols.
  This table is saved in `../data/reference/ensembl_symbol.tsv`.
  - Second, it identifies solid tumor bulk libraries of interest and obtains their `quant.sf` files from S3.
  These files will be used to calculate TPM values as input to deconvolution.
    - Files are saved in `../data/salmon-quant-files/<project id>/<sample id>/<library id>`.
  - To run this script as a member of the Data Lab, you may need to preface with `op run --` if you are using 1Password to manage credentials:
     ```sh
     op run -- Rscript prepare-data-files.R
     ```
- `calculate-tpm.R`: This script creates a TSV of TPM values for all samples in a given project.
Ensembl ids are also converted to gene symbols, where TPM values for duplicate gene symbols are summed.
TSV files are exported to `../data/tpm/<project id>-tpm.tsv`.
