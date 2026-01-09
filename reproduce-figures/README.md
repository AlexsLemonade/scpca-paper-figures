This directory contains a script and certain data files necessary for reproducing manuscript figures, tables, and analyses.
It contains the following contents:

* The script `prepare-scpca-portal-data.R` which prepares data obtained from the ScPCA Portal and AWS S3 needed to reproduce figures, tables, and the bulk RNA-Seq analysis.
Please see instructions below for usage.
  * This script uses functions provided in `utils.R`.
* TSV files, `FigS1B-gene-expression-data.tsv` and `FigS1C_FigS1D-cell-metrics-data.tsv`, which are used as input to the script `../scripts/FigS1B-D_method-metrics-comparison.R` to create Figures S1-D.
  * Note that this script expects these TSV files to be present in this directory; please do not move them.


**CAUTION: You may need up to 300 GB of free storage space for downloaded ScCPA Portal data and associated data preparation.**
As part of data download as described below, you will need to download ZIP files from the ScPCA Portal.
The 300 GB space requirement assumes ZIP files are retained after extraction; deleting ZIP files immediately will reduce the required space.

## Instructions

To reproduce figures and analyses in this repository, you will need to set up your software environment and download ScPCA Portal data.

### Step 1: Set up the `renv` environment

This repository uses [`renv`](https://rstudio.github.io/renv/) to manage R packages and dependencies.
To set up the `renv` environment, open the the `scpca-paper-figures.Rproj` project in RStudio and run `renv::restore()` to install all packages in the provided `renv.lock` file.
Note that this project uses R version `4.4.0`.

Please refer to the [Introduction to `renv` vignette](https://rstudio.github.io/renv/articles/renv.html) for more information about getting started with `renv`.

### Step 2: Download files from the ScPCA Portal

You will need to download all single-cell data (in the `SingleCellExperiment (R)` data format) and metadata from the ScPCA Portal, along with one stand-alone project.
Take the following steps:

1. Navigate to the `Portal-wide downloads` page of the ScPCA Portal: <https://scpca.alexslemonade.org/portal-wide>
2. From this page, download the `Sample Metadata Download`.
  * This will download a ZIP file called `portal-wide_metadata_<accessed-date>.zip`.
    You should unzip this file, which will create a directory `portal-wide_metadata_<accessed-date>`.
    Note the location of this directory for script input.
    You can delete the ZIP file at this point if you prefer.
3. From this page, download the `SingleCellExperiment (R) Download`.
  * Do _not_ click "Merge samples into one project per sample."
  * This will download a ZIP file called `portal-wide_single-cell-experiment_<accessed-date>.zip`.
    You should unzip this file, which will create a directory `portal-wide_single-cell-experiment_<accessed-date>`.
    Note the location of this directory for script input.
    You can delete the ZIP file at this point if you prefer.
4. Next, navigate to the Portal project page for `SCPCP000004`: <https://scpca.alexslemonade.org/projects/SCPCP000004>
  * Click `Download Now` to download this project as a single merged object, selecting the following settings:
    * Modality: `Single-cell`
    * Data Format: `SingleCellExperiment (R)`
    * **Check the box to** `Merge samples into 1 object`
  * This will download a ZIP file called `SCPCP000004_single-cell-experiment_<accessed-date>.zip`.
    You should unzip this file, which will create a directory `SCPCP000004_single-cell-experiment_<accessed-date>`.
    Note the location of this directory for script input.
    You can delete the ZIP file at this point if you prefer.

### Step 3: Run the `prepare-scpca-portal-data.R` script

Now, you can run the `prepare-scpca-portal-data.R` script specifying the following input arguments:

* `--portal_metadata_dir`: Path to the directory containing the downloaded Portal-wide metadata (Step 2 above)
* `--portal_projects_dir`: Path to the directory containing the downloaded Portal-wide single-cell data (Step 3 above)
* `--merged_project_dir`:  Path to the directory containing the downloaded `SCPCP000004` project as a single merged object (Step 4 above)

```sh
Rscript prepare-scpca-portal-data.R \
  --portal_metadata_dir <path to directory with portal-wide metadata > \
  --portal_download_dir <path to directory with portal-wide single-cell data download> \
  --merged_project_dir <path to directory with merged project SCPCP000004>
```

This script will run for roughly 75 minutes, and it will create three directories within this repository with files organized as code expects:

* `../s3_files`: This directory will contain ScPCA data, metadata, and reference files needed to reproduce figures, tables, and the bulk expression analysis
* `../analysis/pseudobulk-bulk-prediction/data/scpca_data`: This directory will contain ScPCA data files needed to reproduce the bulk expression analysis
* `../analysis/pseudobulk-bulk-prediction/data/references`: This directory will contain reference files needed to reproduce the bulk expression analysis

Note that this script will not overwrite existing output directories without the `--overwrite` flag.
For more options, use the `--help` flag.

After running the script, you can remove all downloaded portal ZIP files and directories as desired.

### Step 5: Run repository code

Once `prepare-scpca-portal-data.R` script has successfully completed running, you will be able to run the following repository scripts to reproduce figures and results:

* [`generate-figures-tables.sh`](../generate-figures-tables.sh): This script creates figure panels and tables used in the manuscript by running all R scripts in `../scripts/`.
  * As such you can now run individual scripts in `../scripts/`; however, the scripts in `../scripts/figure_setup/` are intended for internal Data Lab use only.
* [`analysis/pseudobulk-bulk-prediction/run-prediction.sh`](../analysis/pseudobulk-bulk-prediction/run-prediction.sh): This script runs the bulk/pseudobulk expression analysis, which may take several hours to complete.
See the [`analysis` README.md](../analysis/pseudobulk-bulk-prediction/README.md) for more details.