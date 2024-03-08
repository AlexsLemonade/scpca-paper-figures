# scpca-paper-figures

This repo contains the figures and tables included in the ScPCA manuscript.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Summary of figures and tables](#summary-of-figures-and-tables)
- [Generating figures and tables](#generating-figures-and-tables)
- [Sample info](#sample-info)
- [Color palettes](#color-palettes)
- [Renv](#renv)
- [Contributing](#contributing)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Summary of figures and tables

Below is a summary of all figures and tables in the paper.

**Figure 1**

- A. Bar chart summarizing the types of diagnoses found on the Portal.
- B. Bar chart summarizing the types of modalities found on the Portal.
- C. Example project card as shown on the Portal.

**Figure 2**

- A. Overview of main workflow.
- B-G. Simplified versions of plots in the main QC report.

**Figure 3**

- A. Illustration of individual sample download folder.
- B. Illustration of merged project download folder.
- C. Overview of merged workflow.
- D. Example UMAPs found in merged report.

**Figure 4**

- A. Overview of cell type annotation workflow.
- B. Comparison of cell type annotations obtained using `SingleR` and `CellAssign`.

**Supplemental Figure 1**

- A. Comparison of run time and peak memory usage between Alevin-fry and Cell Ranger.
- B. Total UMI/cell between Alevin-fry and Cell Ranger.
- C. Total genes detected/cell between Alevin-fry and Cell Ranger.
- D. Correlation of mean gene expression between Alevin-fry and Cell Ranger.

**Supplemental Figure 2**

- A. Overview of cell type annotation workflow.
- B-D. Simplified versions of plots found in the ADT section of the main QC report.
- E. Overview of multiplexed library workflow.

**Supplemental Figure 3**

- A. Overview of bulk RNA-seq workflow.
- B. Overview of spatial transcriptomics workflow.


**Supplemental Figure 4**

- A. Diagnostic plot for cell type annotations from `SingleR`.
- B. Diagnostic plot for cell type annotations from `CellAssign`.


**Supplemental Figure 5**

Comparison between submitter provided annotations and automated annotations from `SingleR` and `CellAssign`.

**Supplemental Figure 6**

Comparison of delta median statistic obtained from running `SingleR` with different `celldex` references.


**Table S1**

Summary of libraries and types of libraries found on the Portal.

**Table S2**

List of references used for each project on the Portal with `CellAssign`, including the list of organs used to create the reference.


## Generating figures and tables

The `figures` and `tables` folders contain the most up-to-date version of each of the figures and tables, respectively.
The `scripts` directory contains all scripts used to create the figures and tables.
See the [`README` for the `scripts` directory](./scripts/README.md) for more information on figure and table scripts.

To generate all figures and tables, run the script [`generate-figures-tables.sh`](generate-figures-tables.sh) as:

```sh
bash generate-figures-tables.sh
```

Note that this script assumes that the `s3_files` directory has been populated with relevant data files.
These files can be obtained by first running the figure setup scripts, which currently require S3 bucket access as a Data Lab member.
Setup scripts can be run as:

```sh
Rscript scripts/figure_setup/sync-metadata.R
Rscript scripts/figure_setup/sync-data-files.R
Rscript scripts/figure_setup/sync-reference-files.R
```

If you have setup `1Password` to handle your AWS credentials, you will need to prefix those lines with `op run --`:

```sh
op run -- Rscript scripts/figure_setup/sync-metadata.R
op run -- Rscript scripts/figure_setup/sync-data-files.R
op run -- Rscript scripts/figure_setup/sync-reference-files.R
```


## Sample info

The `sample-info` folder contains metadata files used to create figures and tables.

1. `diagnosis-groupings.tsv`: This tsv file contains one row per `submitted_diagnosis` associated with samples on the ScPCA Portal.
For each `submitted_diagnosis`, a `diagnosis_group` is assigned.

2. `disease-timing.tsv`: This tsv file contains one row per `submitted_disease_timing` associated with samples on the ScPCA Portal.
For each `submitted_disease_timing`, a `standardized_disease_timing` is assigned.

3. `project-whiteliest.txt`: This file contains a list of all projects that are currently active on the ScPCA Portal.

## Color palettes

The `palettes` folder contains any palettes used in generating the figures.

1. `diagnosis-group-palette.tsv`: This is the palette used to color the `diagnosis_group` for each sample.
2. `suspension-palette.tsv`: This is the palette used to color libraries by `Single-cell` or `Single-nuclei`.
3. `method-palette.tsv`: This is the palette used to color by quantification method used, either `Alevin-fry` or `Cell Ranger`.

## Manuscript numbers

The `manuscript-numbers` folder contains tables with total sample counts referenced when writing the manuscript.
These tables are not included in the final manuscript and were created using `scripts/Fig1A_sample-disease-barchart.R`.

## Renv

Package dependencies for scripts used in this repo are managed using [`renv`](https://rstudio.github.io/renv/index.html).
For `renv` to work as intended, you'll need to work within the `scpca-paper-figures.Rproj` project in RStudio.
You may need to run `renv::restore()` upon opening the project to ensure the `renv.lock` file is synced with the project library.


## Contributing

When developing new scripts, you may need to install or use new R packages.
Each time you install or use new packages, you will want to run `renv::snapshot()` to update the `renv.lock` file with any added package and dependencies necessary to run the analyses and scripts in this repo.

In addition, this repository uses the [`parsable-r`](https://lorenzwalthert.github.io/precommit/articles/available-hooks.html#parsable-r) pre-commit hook to ensure R scripts are parsable.
To use this hook, first ensure that that the `pre-commit` package is installed on your system; you can install it with your favorite method (`pip install pre-commit` or `conda install pre-commit`).
Then, run `pre-commit install` in the `scpca-paper-figures` directory to enable pre-commit hooks in this repository.
This will install the hooks in the `.git/hooks` directory, and they will be run automatically when you commit changes.
If the hook fails, the commit will be aborted, and you will need to fix the errors and re-commit.
