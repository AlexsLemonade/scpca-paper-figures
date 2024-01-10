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

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Summary of figures and tables

Below is a summary of all figures in the paper.

**Figure 1**

A. Bar chart summarizing the types of diagnoses found on the Portal.
B. Bar chart summarizing the types of modalities found on the Portal.

**Table S1**

Summary of libraries and types of libraries found on the Portal.

## Generating figures and tables

The `figures` and `tables` folders contain the most up-to-date version of each of the figures and tables, respectively.
The `scripts` directory contains all scripts used to create the figures and tables.
See the [`README` for the `scripts` directory](./scripts/README.md) for more information on generating figures.

## Sample info

The `sample_info` folder contains metadata files used to create figures and tables.

1. `diagnosis-groupings.tsv`: This tsv file contains one row per `submitted_diagnosis` associated with samples on the ScPCA Portal.
For each `submitted_diagnosis`, a `diagnosis_group` is assigned.

2. `disease-timing.tsv`: This tsv file contains one row per `submitted_disease_timing` associated with samples on the ScPCA Portal.
For each `submitted_disease_timing`, a `standardized_disease_timing` is assigned.

3. `project-whiteliest.txt`: This file contains a list of all projects that are currently active on the ScPCA Portal.

## Color palettes

The `palettes` folder contains any palettes used in generating the figures.

1. `diagnosis-group-palette.tsv`: This is the palette used to color the `diagnosis_group` for each sample.
2. `suspension-palette.tsv`: This is the palette used to color libraries by `Single-cell` or `Single-nuclei`.

## Renv

Package dependencies for scripts used in this repo are managed using [`renv`](https://rstudio.github.io/renv/index.html).
For `renv` to work as intended, you'll need to work within the `scpca-paper-figures.Rproj` project in RStudio.
You may need to run `renv::restore()` upon opening the project to ensure the `renv.lock` file is synced with the project library.

Each time you install or use new packages, you will want to run `renv::snapshot()` to update the `renv.lock` file with any added package and dependencies necessary to run the analyses and scripts in this repo.
