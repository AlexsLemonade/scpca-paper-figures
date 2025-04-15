This directory contains scripts and results for predicting bulk expression from matched pseudobulk expression.


The analysis takes the following steps:

* Transform bulk count matrices with `DESeq2`
* Prepare pseudobulk matrices and transform with `DESeq2`
* Prepare consensus cell types gene sets for use in overrepresentation analysis
* Construct linear models predicting bulk from pseudobulk expression
* Perform overrepresentation analysis on residual outliers

## Running the analysis

To run the analysis, you will first need to prepare input data files as described below.
Then, you can run the analysis bash script:

```sh
bash run-prediction.sh
```

## Instructions to prepare input data

To prepare and obtain all files needed for analysis input, please follow the instructions provided in [`reproduce-figures/README.md`](../../reproduce-figures/README.md).

### Instructions for Data Lab members

To prepare data for analysis, you will need to take the following steps.
Note that you may need to preface these lines with `op run --` if you are using 1Password to manage AWS credentials.

1. Sync files to `./data/`:

```sh
Rscript scripts/sync-data-files.R \
  --output_dir "data/scpca_data" \
  --reference_dir "data/references" \
  --map_file "data/bulk-library-sample-ids.tsv"
```

1. Ensure additional input files have been synced to `../../s3_files`:

```sh
# These should be run from top-level of the repository
Rscript scripts/figure_setup/sync-consensus-celltype-results.R
Rscript scripts/figure_setup/sync-metadata.R
```
