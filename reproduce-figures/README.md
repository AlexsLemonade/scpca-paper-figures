This directory contains a script and certain data files necessary for reproducing manuscript figures, tables, and analyses.
It contains the following contents:

* The script `prepare-scpca-portal-data.R` prepares data obtained from the ScPCA Portal and AWS S3 needed to reproduce figures, tables, and the bulk RNA-Seq analysis.
Please see instructions below for usage.
  * This script uses functions provided in `utils.R`
  * This script also uses the `./scratch` directory during file processing
* TSV files `FigS1B-gene-expression-data.tsv` and `FigS1C_FigS1D-cell-metrics-data.tsv` are used as input to the script `../scripts/FigS1B-D_method-metrics-comparison.R` which creates Figures S1-D.
  * Note that this script expects these TSV files to be present in this directory; please do not move them.

## Instructions

To reproduce figures and analyses in this repository, you will need to organize and prepare all data needed for input.
Follow the steps below to do so:

### Step 1: Set up the `renv` environment

This repository uses [`renv`](https://rstudio.github.io/renv/) to manage R packages and dependencies.
To set up the `renv` environment, open the the `scpca-paper-figures.Rproj` project in RStudio and run `renv::restore()` to install all packages in the provided `renv.lock` file.
Note that this project uses R version `4.4.0`.

Please refer to the [Introduction to `renv` vignette](https://rstudio.github.io/renv/articles/renv.html) for more information about getting started with `renv`.

### Step 2: Download files from the ScPCA Portal

You will need to download the following files from the ScPCA Portal <https://scpca.alexslemonade.org/>:

1. All ScPCA projects listed in the [Project Whitelist file](../sample-info/project-whitelist.txt)

  * On each project page (<https://scpca.alexslemonade.org/projects/SCPCPXXXXXX>), download each project with the following options selected:
    * Ensure the selected modality is `Single-cell`
    * Ensure the selected Data Format is `SingleCellExperiment (r)`
    * Do _not_ click to merge all objects into the same sample
    * If the project contains multiplexed samples, you can exclude those samples (but it will not affect the code if they are included)
  * These projects will download as zip files; you should leave them compressed and named as is.

2. [TBD] The merged project for `SCPCP000003` only, again in `SingleCellExperiment (r)` format.

3. The portal-wide metadata TSV file.
You can download this file using the "Get All Sample Metadata" button on the top-right of the portal homepage.


### Step 3: Organize your downloaded files

Next, you will organize this data for input to the `prepare-scpca-portal-data.R` script:

* Place all project ZIP files into a single directory
* Optionally, you can also store the merged `SCPCP000003` project ZIP file and the portal metadata TSV in this directory, but it is not necessary.

### Step 4: Run the `prepare-scpca-portal-data.R` script

Now, you can run the `prepare-scpca-portal-data.R` script as follows:

**Caution: Running this script requires at least 170 GB free storage space and will take roughly 90 minutes.**

```sh
  Rscript prepare-scpca-portal-data.R \
    --portal_projects_dir <path to directory with all project zip files> \
    --merged_sce_path <path to merged SCE file> \
    --portal_metadata_path <path to portal-wide metadata TSV>
```

This script will create three directories within this repository with files organized as code expects:

* `../s3_files`: This directory will contain ScPCA data, metadata, and reference files needed to reproduce figures, tables, and the bulk expression analysis
* `../analysis/pseudobulk-bulk-prediction/data/scpca_data`: This directory will contain ScPCA data files needed to reproduce the bulk expression analysis
* `../analysis/pseudobulk-bulk-prediction/data/references`: This directory will contain reference files needed to reproduce the bulk expression analysis

Note that this script will not overwrite existing output directories without the `--overwrite` flag.

### Step 5: Run repository code

Once `prepare-scpca-portal-data.R` script has successfully completed running, you will be able to run the following repository scripts to reproduce figures and results

* [`generate-figures-tables.sh`](../generate-figures-tables.sh): This script creates figure panels and tables used in the manuscript by running all R scripts in `scripts/`.
  * As such you can now run individual scripts in `scripts/`; however, the scripts in `scripts/figure_setup/` are intended for internal Data Lab use only.
* [`analysis/pseudobulk-bulk-prediction/run-prediction.sh`](../analysis/pseudobulk-bulk-prediction/run-prediction.sh): This script runs the bulk/pseudobulk expression analysis, which may take several hours to complete.