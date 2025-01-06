This directory contains scripts used in the bulk deconvolution analysis.

- `sync-salmon-output.R`: This script identifies solid tumor bulk libraries of interest and obtains their `quant.sf` files from S3.
These files will be used to calculate TPM values as input to deconvolution.
This script saves all files in `../data/salmon-quant-files/<project id>/<sample id>/<library id>`; note that the `data` directory is ignored by Git.
 - To run this script as a member of the Data Lab, you may need to preface with `op run --` if you are using 1Password to manage credentials:
   ```sh
   op run -- Rscript sync-salmon-output.R
   ```
