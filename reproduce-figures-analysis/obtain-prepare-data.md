# Preparing input files for figures and analysis

This file contains instructions for obtaining and preparing data files to reproduce figures, tables, and analyses in `scpca-paper-figures`.

In addition to files provided already in the repository, you will also need the following:

* Files from the ScPCA Portal
  * All single-cell files from the portal shoud be downloaded in `SingleCellExperiment (R)` format
* Files from the public AWS S3 bucket `s3://scpca-references`
  * To obtain files from S3, you will need to use [the `awscli` tool](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Data for figure and table scripts](#data-for-figure-and-table-scripts)
  - [File organization](#file-organization)
  - [Data files from the ScPCA Portal](#data-files-from-the-scpca-portal)
  - [Reference files](#reference-files)
  - [ScPCA metadata files](#scpca-metadata-files)
  - [TODO: Files to reproduce Figures S1B-D](#todo-files-to-reproduce-figures-s1b-d)
  - [TODO: Consensus cell type files](#todo-consensus-cell-type-files)
- [Bulk and pseudobulk RNASeq analysis](#bulk-and-pseudobulk-rnaseq-analysis)
  - [Obtaining files in `references`](#obtaining-files-in-references)
  - [Obtaining files in `scpca_data`](#obtaining-files-in-scpca_data)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->



## Data for figure and table scripts

This section provides instructions on how to obtain and prepare data files to run code in the `figure_scripts/` directory.

### File organization

To faciliate reproducing figures, we recommend organizing files into this structure in a directory called `s3_files`, to be stored at the top-level of the repository.
Instructions on how to obtain each of these files are given in the following sections.

```console
├── SCPCP000003
│   └── SCPCP000003_merged.rds
├── SCPCS000001
│   ├── SCPCL000001_filtered.rds
│   ├── SCPCL000001_processed.rds
│   └── SCPCL000001_unfiltered.rds
├── SCPCS000216
│   ├── SCPCL000290_filtered.rds
│   ├── SCPCL000290_processed.rds
│   └── SCPCL000290_unfiltered.rds
├── SCPCS000251
│   └── SCPCL000498_processed.rds
├── SCPCS000264
│   └── SCPCL000490_processed.rds
├── cell-type-consensus-results
│   └──....
├── celltype_results
│   ├── SCPCL000001_processed.rds
│   ├── SCPCL000002_processed.rds
│   ├── SCPCL000004_processed.rds
│   ├── SCPCL000296_processed.rds
│   ├── SCPCL000297_processed.rds
│   ├── SCPCL000298_processed.rds
│   ├── SCPCL000478_processed.rds
│   ├── SCPCL000484_processed.rds
│   └── SCPCL000488_processed.rds
├── reference_files
│   ├── Homo_sapiens.GRCh38.104.mitogenes.txt
│   ├── Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv
│   └── singler_models
│       ├── BlueprintEncodeData_celldex_1-10-1_model.rds
│       ├── DatabaseImmuneCellExpressionData_celldex_1-10-1_model.rds
│       ├── HumanPrimaryCellAtlasData_celldex_1-10-1_model.rds
│       └── MonacoImmuneData_celldex_1-10-1_model.rds
├── scpca-library-metadata.tsv
├── scpca-project-celltype-metadata.tsv
└── scpca-sample-metadata.tsv
```

### Data files from the ScPCA Portal

The table below provides an overview of RDS files with `SingleCellExperiment` that need to be downloaded from the ScPCA Portal to reproduce figures:

| File                                                | Link to ScPCA Portal Project                           | Notes                                                                         |
| --------------------------------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------------------- |
| `SCPCL000001_<unfiltered, filtered, processed>.rds` | <https://scpca.alexslemonade.org/projects/SCPCP000001> | Download files for sample `SCPCS000001`                                       |
| `SCPCL000002_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000001> | Download files for sample `SCPCS000002`                                       |
| `SCPCL000004_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000001> | Download files for sample `SCPCS000004`                                       |
| `SCPCP000003_merged.rds`                            | <https://scpca.alexslemonade.org/projects/SCPCP000003> | Download the full project, selecting the option `Merge samples into 1 object` |
| `SCPCL000478_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000005> | Download files for sample `SCPCS000252`                                       |
| `SCPCL000484_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000005> | Download files for sample `SCPCS000258`                                       |
| `SCPCL000488_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000005> | Download files for sample `SCPCS000262`                                       |
| `SCPCL000290_<unfiltered, filtered, processed>.rds` | <https://scpca.alexslemonade.org/projects/SCPCP000007> | Download files for sample `SCPCS000216`                                       |
| `SCPCL000296_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000007> | Download files for sample `SCPCS000222`                                       |
| `SCPCL000297_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000007> | Download files for sample `SCPCS000223`                                       |
| `SCPCL000298_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000007> | Download files for sample `SCPCS000224`                                       |
| `SCPCL000298_processed.rds`                         | <https://scpca.alexslemonade.org/projects/SCPCP000007> | Download files for sample `SCPCS000224`                                       |



### Reference files

All files in the `reference_files` directory are available from the public AWS S3 bucket `s3://scpca-references`.
They can be obtained using the `awscli` tool with the following commands:

```sh
aws s3 cp s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv . --no-sign-request
aws s3 cp s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt . --no-sign-request
aws s3 cp s3://scpca-references/celltype/singler_models/ . --recursive --exclude "*" --include "*.rds" --no-sign-request
```


### ScPCA metadata files

There are three metadata files you will need:

* `scpca-project-celltype-metadata.tsv`
* `scpca-library-metadata.tsv`
* `scpca-sample-metadata.tsv`

A version of `scpca-project-celltype-metadata.tsv` which you can use to reproduce figures is provided in TODO.

To obtain the other two metadata files, we offer a helper script in TODO to create versions of these files which you can use to reproduce figures.
Follow these instructions to create these two metadata files:

* From the [ScPCA Portal](https://scpca.alexslemonade.org/), download the full portal metadata using the `Get All Sample Metadata` button on the top-right of the page
* Run the helper script as follows:

```sh
# TODO!!
Rscript prepare-metadata-files.tsv --all_sample_metadata <path to that downloaded file>
```

This will create files named `scpca-library-metadata.tsv` and `scpca-sample-metadata.tsv` which you can use to recreate figures.
These files will be created in the same directory from which you run the script.


### TODO: Files to reproduce Figures S1B-D

Forthcoming section: Here, we can describe the benchmarking TSV files (https://github.com/AlexsLemonade/scpca-paper-figures/issues/216)

### TODO: Consensus cell type files

Forthcoming section: Here, we can describe obtaining consensus cell type files which is TBD.

## Bulk and pseudobulk RNASeq analysis

Code for this analysis is provided in the directory `analysis/pseudobulk-bulk-prediction`.
From this directory, you can use the `run-prediction.sh` script to run the analysis and generate results which are used to as input to scripts in `figure_scripts` that create Figure 6 and Figure S8.

There are several files you will need to obtain from both the ScPCA Portal and the public AWS S3 bucket to enable running the analysis.
Data should be stored in `analysis/pseudobulk-bulk-prediction/data` according to this organization:

```console
├── references
│   ├── bone-and-soft-tissue_PanglaoDB_2020-03-27.tsv
│   ├── brain-compartment_PanglaoDB_2020-03-27.tsv
│   └── kidney-compartment_PanglaoDB_2020-03-27.tsv
└── scpca_data
    ├── SCPCP000001
    │   ├── SCPCP000001_bulk_quant.tsv
    │   ├── SCPCSXXXXXX
    │   │   └── SCPCLXXXXXX_processed.rds
    │   └── ...
    ├── SCPCP000002
    │   ├── SCPCP000002_bulk_quant.tsv
    │   ├── SCPCSXXXXXX
    │   │   └── SCPCLXXXXXX_processed.rds
    │   └── ...
    ├── SCPCP000006
    │   ├── SCPCP000006_bulk_quant.tsv
    │   ├── SCPCSXXXXXX
    │   │   └── SCPCLXXXXXX_processed.rds
    │   └── ...
    ├── SCPCP000009
    │   ├── SCPCP000009_bulk_quant.tsv
    │   ├── SCPCSXXXXXX
    │   │   └── SCPCLXXXXXX_processed.rds
    │   └── ...
    └── SCPCP000017
        ├── SCPCP000017_bulk_quant.tsv
        ├── SCPCSXXXXXX
        │   └── SCPCLXXXXXX_processed.rds
        └── ...
```

### Obtaining files in `references`

Files in the `references` directory can be obtained from the public `scpca-references` AWS S3 bucket using [the `awscli` tool](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

Use the following commands to download files into the current directory:

```sh
aws s3 cp s3://scpca-references/celltype/cellassign_references/bone-and-soft-tissue_PanglaoDB_2020-03-27.tsv . --no-sign-request
aws s3 cp s3://scpca-references/celltype/cellassign_references/brain-compartment_PanglaoDB_2020-03-27.tsv . --no-sign-request
aws s3 cp s3://scpca-references/celltype/cellassign_references/kidney-compartment_PanglaoDB_2020-03-27.tsv . --no-sign-request
```

### Obtaining files in `scpca_data`

All ScPCA data files come from the ScPCA Portal projects `SCPCP000001`, `SCPCP000002`, `SCPCP000006`, `SCPCP000009`, and `SCPCP000017`.
To obtain the `_processed.rds` files associated with samples of interest, we recommend taking the following steps:

* Navigate to the ScPCA Portal project page of interest, e.g. <https://scpca.alexslemonade.org/projects/SCPCP000001> for project `SCPCP000001`.
* Click the "Download Project" button to download all samples for the project.
Do _not_ click "Merge all samples into 1 object", and be sure to use the `SingleCellExperiment (R)` format when downloading.
* This will download both project's `SingleCellExperiment` files, organized by sample, as well as the project's bulk raw counts matrix
* Once downloaded, you will then need to _remove_ directories for single-cell samples which are not used in this analysis.
    * We provide a list in `bulk-analysis-samples.tsv` of all sample ids which _are used_; other samples not listed here shoud be removed.
    * Note that project `SCPCP000009` contains samples processed with multiplexed libraries, whose directory names include underscores; none of these are included in analysis, and all such directories should be removed.
