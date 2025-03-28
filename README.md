# scpca-paper-figures

This repo contains the figures and tables included in the ScPCA manuscript.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Summary of figures and tables](#summary-of-figures-and-tables)
- [Generating figures and tables](#generating-figures-and-tables)
  - [Instructions to prepare data](#instructions-to-prepare-data)
  - [Instructons for Data Lab members](#instructons-for-data-lab-members)
- [Additional repository contents](#additional-repository-contents)
  - [Sample info](#sample-info)
  - [Color palettes](#color-palettes)
  - [Manuscript numbers](#manuscript-numbers)
  - [Nextflow logs](#nextflow-logs)
  - [Analysis](#analysis)
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

**Figure 5**

- A. Dot plot summarizing marker gene expression in consensus cell types in Brain and CNS samples.
- B. Bar plot summarizing the percentage of cells annotated as each consensus cell type in High-grade and Low-grade glioma samples.
- C. Bar plot summarizing the percentage of cells annotated with specific immune cell types in High-grade and Low-grade glioma samples.

**Figure 6**

- A. Scatterplots of the relationship between bulk and pseudobulk counts for relevant Brain and CNS projects.
- B. Bar plots of odds ratios from overrepresentation analysis of bulk expression data for relevant Brain and CNS projects.

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

Comparison of delta median statistic obtained from running `SingleR` with different `celldex` references.

**Supplemental Figure 5**

- A. UMAP displaying cell type annotations from `CellAssign` for an example ScPCA library.
- B. Heatmap comparing submitter provided annotations to `CellAssign` and `SingleR` annotations for an example ScPCA library.

**Supplemental Figure 6**

- A. Dot plot summarizing marker gene expression in consensus cell types in Leukemia samples.
- B. Dot plot summarizing marker gene expression in consensus cell types in Sarcoma samples.
- C. Dot plot summarizing marker gene expression in consensus cell types in Other solid tumor samples.

**Supplemental Figure 7**

- A. Bar plot summarizing the percentage of cells annotated as each consensus cell type in Leukemia samples.
- B. Bar plot summarizing the percentage of cells annotated as each consensus cell type in Sarcoma samples.
- C. Bar plot summarizing the percentage of cells annotated as each consensus cell type in Other solid tumor samples.

**Supplemental Figure 8**

- A. Scatterplots of the relationship between bulk and pseudobulk counts for projects not shown in Figure 6A.
- B. Bar plots of odds ratios from overrepresentation analysis of bulk expression data for projects not shown in Figure 6B.


**Table S1**

Summary of libraries and types of libraries found on the Portal.

**Table S2**

List of references used for each project on the Portal with `CellAssign`, including the list of organs used to create the reference.


## Generating figures and tables

The `figures/` and `tables/` folders contain the most up-to-date version of each of the figures and tables, respectively.
The `scripts/` folder contains all scripts used to create the figures and tables.
See the [`README` for the `scripts` folder](./scripts/README.md) for more information on figure and table scripts.


The [`generate-figures-tables.sh`](generate-figures-tables.sh) script can be used to prepare all figures and tables.
The script assumes that a folder called `s3_files` has been populated with relevant data files.
Therefore, **before running this script** you will first need to prepare these input data files, as described in the sections below.

Then, run the figure generation script as:

```sh
bash generate-figures-tables.sh
```

### Instructions to prepare data

If you are not a member of the Data Lab, please follow the instructions provided in `reproduce-figures-analysis/obtain-prepare-data.md` to obtain and prepare additional data needed to regenerate figures and tables.

### Instructons for Data Lab members

To prepare data for figure and table generation, will need to run the figure setup scripts:

```sh
Rscript scripts/figure_setup/sync-metadata.R
Rscript scripts/figure_setup/sync-data-files.R
Rscript scripts/figure_setup/sync-reference-files.R
Rscript scripts/figure_setup/assign-brain-classifications.R
Rscript scripts/figure_setup/prepare-sample-whitelist.R
```

If you have setup `1Password` to handle your AWS credentials, you will need to prefix scripts beginning with `sync-` with `op run --`, specifically:

```sh
op run -- Rscript scripts/figure_setup/sync-metadata.R
op run -- Rscript scripts/figure_setup/sync-data-files.R
op run -- Rscript scripts/figure_setup/sync-reference-files.R
```

## Additional repository contents

### Sample info

The `sample-info/` folder contains metadata files used to create figures and tables.

* `brain-classifications-no-multiplexed.tsv`: This tsv file classifies brain-related diagnoses in the ScPCA Portal into "High-grade glioma" and "Low-grade glioma" for plotting.
* `celltype-reference-metadata.tsv`: This tsv file contains information about references used for CellAssign cell type annotation on the ScPCA Portal.
* `diagnosis-groupings.tsv`: This tsv file contains one row per `submitted_diagnosis` associated with samples on the ScPCA Portal.
For each `submitted_diagnosis`, a `diagnosis_group` is assigned.
* `disease-timing.tsv`: This tsv file contains one row per `submitted_disease_timing` associated with samples on the ScPCA Portal.
For each `submitted_disease_timing`, a `standardized_disease_timing` is assigned.
* `project-whitelist.txt`: This file contains a list of all projects that are currently active on the ScPCA Portal.
* `sample-whitelist.txt`: This file contains a list of all samples that are currently active on the ScPCA Portal.
* `scpca-project-celltype-metadata.tsv`: This tsv file provides the specific CellAssign cell typing reference used for each ScPCA project.


### Color palettes

The `palettes/` folder contains any palettes used in generating the figures.

* `diagnosis-group-palette.tsv`: This is the palette used to color the `diagnosis_group` for each sample.
* `disease-timing-palette.tsv`: This is the palette used to color the `disease_timing` for each sample.
* `immune-palette.tsv`: This is the palette used to color certain immune cell types from the overall consensus cell types.
* `suspension-palette.tsv`: This is the palette used to color libraries by `Single-cell` or `Single-nuclei`.
* `method-palette.tsv`: This is the palette used to color by quantification method used, either `Alevin-fry` or `Cell Ranger`.
* `validation-group-palette.tsv`: This is the palette used to color validation groups used to assess consensus cell types.

### Manuscript numbers

The `manuscript-numbers` folder contains tables with total sample counts referenced when writing the manuscript.
These tables are not included in the final manuscript and were created as follows:

* `bulk-analysis-counts.tsv` was created by `../scripts/Fig6-FigS8_bulk-analysis.R`
* `diagnosis-group-counts.tsv` and `disease-timing-counts.tsv` were created by `../scripts/Fig1A_sample-disease-barchart.R`


### Nextflow logs

The `nextflow_logs` folder contains text files with Nextflow log information from running the `scpca-nf` pipeline.
These files are specifically used by `scripts/FigS1A_memory-time-comparison.R` to create Figure 1A.

### Analysis

This `analysis` folder contains analyses comparing bulk and pseudobulk RNA-Seq data.
Please refer to `analysis/README.md` for additional details.

## Renv

Package dependencies for scripts used in this repo are managed using [`renv`](https://rstudio.github.io/renv/index.html).
For `renv` to work as intended, you'll need to work within the `scpca-paper-figures.Rproj` project in RStudio.
You may need to run `renv::restore()` upon opening the project to ensure the `renv.lock` file is synced with the project library.


## Contributing

When developing new scripts, you may need to install or use new R packages.
Each time you install or use new packages, you will want to run `renv::snapshot()` to update the `renv.lock` file with any added package and dependencies necessary to run the analyses and scripts in this repo.

In addition, this repository uses the [`parsable-r`](https://lorenzwalthert.github.io/precommit/articles/available-hooks.html#parsable-r) pre-commit hook to ensure R scripts are parsable.
To use this hook, first ensure that that the `pre-commit` package is installed on your system; you can install it with your favorite method (`pip install pre-commit` or `conda install pre-commit`).
Then, run `pre-commit install` in the `scpca-paper-figures` folder to enable pre-commit hooks in this repository.
This will install the hooks in the `.git/hooks` folder, and they will be run automatically when you commit changes.
If the hook fails, the commit will be aborted, and you will need to fix the errors and re-commit.
