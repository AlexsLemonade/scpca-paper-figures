This directory contains a script and certain data files necessary for reproducing manuscript figures, tables, and analyses.
It contains the following contents:

* The script `prepare-scpca-portal-data.R` prepares data obtained from the ScPCA Portal and AWS S3 needed to reproduce figures, tables, and the bulk RNA-Seq analysis.
Please see instructions below for usage.
  * **CAUTION: Running this script requires at least 170 GB of free storage space**
  * This script uses functions provided in `utils.R`
  * This script also uses the `./scratch` directory during file processing
* TSV files, `FigS1B-gene-expression-data.tsv` and `FigS1C_FigS1D-cell-metrics-data.tsv`, which are used as input to the script `../scripts/FigS1B-D_method-metrics-comparison.R` to create Figures S1-D.
  * Note that this script expects these TSV files to be present in this directory; please do not move them


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

1. All ScPCA projects listed in the [project whitelist file](../sample-info/project-whitelist.txt)
  * Use the project identifiers listed in the project whitelist file to navigate to each project page.
For example, to navigate to `SCPCP000001`, use the URL: <https://scpca.alexslemonade.org/projects/SCPCP000001>.
* Download each project using the `Download Project` button with the following options selected:
    * Ensure the selected modality is `Single-cell`
    * Ensure the selected Data Format is `SingleCellExperiment (r)`
    * Do _not_ click to merge all objects into the same sample
    * If the project contains multiplexed samples, you can exclude those samples (but it will not affect the code if they are included)
  * These projects will download as zip files; you should leave them in this compressed format, and _do not rename them_.

2. The portal-wide metadata TSV file
  * You can download this single TSV file using the "Get All Sample Metadata" button on the top-right of the portal homepage.


### Step 3: Organize your downloaded files

Place all downloaded ScPCA project files in the same directory.
Optionally, you can also store the portal metadata TSV in this directory, but it is not necessary.

### Step 4: Run the `prepare-scpca-portal-data.R` script

**CAUTION: Running this script requires at least 170 GB of free storage space.**

Now, you can run the `prepare-scpca-portal-data.R` script specifying the following input arguments:

* `--portal_projects_dir`: Path to the directory containing the ScPCA project zip files
* `--portal_metadata_path`: Path to the portal-wide metadata TSV file

```sh
  Rscript prepare-scpca-portal-data.R \
    --portal_projects_dir <path to directory with all ScPCA project zip files> \
    --portal_metadata_path <path to portal-wide metadata TSV>
```

This script will run for roughly 90 minutes, and it will create three directories within this repository with files organized as code expects:

* `../s3_files`: This directory will contain ScPCA data, metadata, and reference files needed to reproduce figures, tables, and the bulk expression analysis
* `../analysis/pseudobulk-bulk-prediction/data/scpca_data`: This directory will contain ScPCA data files needed to reproduce the bulk expression analysis
* `../analysis/pseudobulk-bulk-prediction/data/references`: This directory will contain reference files needed to reproduce the bulk expression analysis

Note that this script will not overwrite existing output directories without the `--overwrite` flag.
For more options, use the `--help` flag.

### Step 5: Run repository code

Once `prepare-scpca-portal-data.R` script has successfully completed running, you will be able to run the following repository scripts to reproduce figures and results:

* [`generate-figures-tables.sh`](../generate-figures-tables.sh): This script creates figure panels and tables used in the manuscript by running all R scripts in `scripts/`.
  * As such you can now run individual scripts in `scripts/`; however, the scripts in `scripts/figure_setup/` are intended for internal Data Lab use only.
* [`analysis/pseudobulk-bulk-prediction/run-prediction.sh`](../analysis/pseudobulk-bulk-prediction/run-prediction.sh): This script runs the bulk/pseudobulk expression analysis, which may take several hours to complete.