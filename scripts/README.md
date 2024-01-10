# Scripts for generating figures and tables

This folder contains all scripts used to generate figures and tables.

## Figure set up

The `figure_setup` folder contains scripts that are required to be run _prior to_ generating figures.

1. `sync-metadata.R`: This script is used to sync any metadata files found on S3 to a local folder.
In particular, the `scpca-sample-metadata.tsv` and `scpca-library-metadata.tsv` will be stored to a folder within the root directory of this repo named `s3_files`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
The `s3_files` folder is ignored in the `.gitignore`, so the first time you go to generate figures you will need to run this script as these files are not available in the repo.
To run the script use the following command:

```sh
op run -- Rscript sync-metadata.R
```

## Generating figures and tables

The following scripts can be used to generate figures:

1. `fig_1a-sample_disease_barchart.R`: This script is used to generate Figure 1A, which includes a summary of the types of diagnoses found in the ScPCA portal.
Before running this script, you must run `figure_setup/sync-metadata.R`.

2. `fig_1b-modality_barchart.R`: This script is used to generate Figure 1B, which includes a summary of the types of modalities found in the ScPCA portal.
Before running this script, you must run `figure_setup/sync-metadata.R`.

3. `table_s1-modality_summary.R`: This script is used to generate supplemental Table 1, which contains a summary of the types of libraries found in each project.
Before running this script, you must run `figure_setup/sync-metadata.R`.
