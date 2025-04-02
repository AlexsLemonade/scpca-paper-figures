# Preparing input files for figures and analysis

This file contains instructions for obtaining and preparing data files to reproduce figures and tables and analyses in `scpca-paper-figures`.

In addition to files provided already in the repository, you will also need the following:

* Files from the ScPCA Portal
  * All single-cell files from the portal shoud be downloaded in `SingleCellExperiment (R)` format
* Files from the public AWS S3 bucket `s3://scpca-references`
  * To obtain files from S3, you will need to use [the `awscli` tool](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [File organization](#file-organization)
- [Data files from the ScPCA Portal](#data-files-from-the-scpca-portal)
- [Reference files](#reference-files)
- [ScPCA metadata files](#scpca-metadata-files)
- [TODO: Files to reproduce Figures S1B-D](#todo-files-to-reproduce-figures-s1b-d)
- [TODO: Consensus cell type files](#todo-consensus-cell-type-files)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


This section provides instructions on how to obtain and prepare data files to run code in the `scripts/` directory.
By preparing data as described in this file, you will not need to run any scripts in `scripts/figure_setup/` which are meant for internal Data Lab use.

## File organization

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
└── scpca-sample-metadata.tsv
```

## Data files from the ScPCA Portal

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



## Reference files

All files in the `reference_files` directory are available from the public AWS S3 bucket `s3://scpca-references`.
They can be obtained using the `awscli` tool with the following commands:

```sh
aws s3 cp s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv . --no-sign-request
aws s3 cp s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt . --no-sign-request
aws s3 cp s3://scpca-references/celltype/singler_models/ . --recursive --exclude "*" --include "*.rds" --no-sign-request
```

## ScPCA metadata files

There are two metadata files you will need, `scpca-library-metadata.tsv` and `scpca-sample-metadata.tsv`.
To create them, we offer a helper script in `prepare-metadata-files.R` to create versions of these files which you can use to reproduce figures.
Follow these instructions to create these two metadata files:

* From the [ScPCA Portal](https://scpca.alexslemonade.org/), download the full portal metadata using the `Get All Sample Metadata` button on the top-right of the page
* Run the helper script as follows:

```sh
# By default, files will be exported a directory called metadata-files/
Rscript prepare-metadata-files.R --all_sample_metadata <path to that downloaded file>

# You can customize this directory with  --output_dir
Rscript prepare-metadata-files.R --all_sample_metadata <path to that downloaded file> --output_dir <output directory>

# If metadata files already exist in the output directory, you will need the --overwrite flag to overwrite them:
Rscript prepare-metadata-files.R --all_sample_metadata <path to that downloaded file> --overwrite
```

This will create files named `scpca-library-metadata.tsv` and `scpca-sample-metadata.tsv` which you can use to recreate figures.
These files will be created in the same directory from which you run the script.


## TODO: Files to reproduce Figures S1B-D

Forthcoming section: Here, we can describe the benchmarking TSV files (https://github.com/AlexsLemonade/scpca-paper-figures/issues/216)

## TODO: Consensus cell type files

Forthcoming section: Here, we can describe obtaining consensus cell type files which is TBD.
