# Scripts for generating figures and tables

This folder contains all scripts used to generate figures and tables.
Before running these scripts, you will need to prepare data for script input:

* If you are a Data Lab staff member, you will need to run scripts described in the **Figure set up** section at the bottom of this document
* You are not a Data Lab staff member, you will need to follow instructions provided in the file `../reproduce-figures-analysis/obtain-prepare-data.md`

## Generating figures and tables

The following scripts can be used to generate figures and tables:

* `Fig1A_sample-disease-barchart.R`: This script is used to generate Figure 1A, which includes a summary of the types of diagnoses found in the ScPCA portal.

* `Fig1B_modality-barchart.R`: This script is used to generate Figure 1B, which includes a summary of the types of modalities found in the ScPCA portal.

* `Fig2B-F_qc-plots.R`: This script is used to generate panels Figure 2B-F, which includes simplified and miniature versions of the plots found in the main QC report included with each sample download.

* `Fig3D_merged-umaps.R`: This script is used to generate Figure 3D, which includes an example of the UMAPs shown in the merged report.

* `Fig4B_singler-cellassign-heatmap.R`: This script is used to generate Figure 4B, which includes an example of the heatmap comparing cell type annotations from `SingleR` to `CellAssign`.

* `Fig5A-S6_consensus-cell-type-dotplots.R`: This script is used to generate Figure 5A and supplemental figure 6A-C.
Figure 5A includes a dot plot summarizing marker gene expression across all consensus cell types in the Brain and CNS samples.
Supplemental Figure 6A-C include dot plots summarizing marker gene expression across all consensus cell types in the other diagnosis groups, leukemia, sarcoma, and other solid tumors.

* `Fig5B-C_brain-barplots.R`: This script is used to generate Figures 5B and 5C.
Figure 5B shows the percentage of each library that is annotated as each of the consensus cell types for all libraries in High-grade and Low-grade glioma samples.
Figure 5C shows the percentage of each library that is annotated as each of the _immune_ consensus cell types, emphasizing myeloid and leukocyte cell types, for all libraries in High-grade and Low-grade glioma samples.

* `Fig6-FigS8_bulk-analysis.R `: This script is used to generate Figures 6 and S8.
Figure 6 shows (A) scatterplots of the relationship between bulk and pseudobulk counts for relevant Brain and CNS projects, and (B) bar plots of odds ratios from overrepresentation analysis of bulk expression data for relevant Brain and CNS projects.
Figure S8 shows (A) scatterplots of the relationship between bulk and pseudobulk counts for projects not shown in Figure 6A, and (B) bar plots of odds ratios from overrepresentation analysis of bulk expression data for projects not shown in Figure 6B.

* `FigS1A_memory-time-comparison.R`: This script is used to generate supplemental Figure 1A, which shows a comparison of total run time and peak memory usage for Cell Ranger and Alevin-fry.
This script uses the trace files found in `nextflow_logs/`.

* `FigS1B-D_method-metrics-comparison.R`: This script is used to generate supplemental figures 1B-D, which compares cell and gene level metrics between libraries quantified using Cell Ranger and Alevin-fry.

* `FigS2B-D_adt-plots.R`: This script is used to generate supplemental Figure 2B-D, which includes simplified and miniature versions of the plots in the ADT section of the main QC report included with each sample download.

* `FigS4_celldex-ref-comparison.R`: This script is used to generate supplemental Figure 4, which compares the delta median statistic calculated from running `SingleR` on a subset of ScPCA libraries with 4 different `celldex` references.

* `FigS5A_cellassign-justification.R`: This script is used to generate supplemental Figure 5A, which includes a UMAP summarizing cell type annotations using `CellAssign`.

* `FigS5B_submitter-celltypes-heatmap.R`: This script is used to generate supplemental Figure 5B, which includes a heatmap comparing the `CellAssign` and `SingleR` annotations to submitter-provided annotations for an example library.

* `FigS7_consensus-bar-plots.R`: This script is used to generate Supplemental figure 7A-C, which show the percentage of each library that is annotated as each of the consensus cell types for all libraries in each diagnosis group, excluding the Brain and CNS samples.

* `TableS1_modality-summary.R`: This script is used to generate supplemental Table 1, which contains a summary of the types of libraries found in each project.

* `TableS2_cellassign-ref-summary.R`: This script is used to generate supplemental Table 2, which includes one row for each ScPCA project on the Portal and the associated diagnoses, reference used from `PanglaoDB`, and list of organs used to construct the reference.

## Old figures

The `old_figures` folder contains the scripts used to generate figures present in previous versions of the manuscript.

* `FigS5B_cellassign-heatmap.R`: This script is used to generate a previous version of supplemental Figure 5B, which included an example heatmap comparing submitter provided annotations to automated annotations from `CellAssign`.

* `FigS6A-B_celltype-diagnostic-plots.R`: This script is used to generate a previous version of supplemental Figure 6A-B, which included examples of the diagnostic plots for annotating cell types with `SingleR` and `CellAssign`.


## Figure set up

The `figure_setup` folder contains scripts to prepare figure data for internal use by the Data Lab.
These scripts are only expected to be run by Data Lab staff.

* `sync-metadata.R`: This script is used to sync any metadata files found on S3 to a local folder.
In particular, the `scpca-sample-metadata.tsv` and `scpca-library-metadata.tsv` will be stored to a folder within the root directory of this repo named `s3_files`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
The `s3_files` folder is included in the `.gitignore`, so the first time you go to generate figures you will need to run this script as these files are not available in the repo.
To run the script use the following command:

```sh
op run -- Rscript sync-metadata.R
```

* `sync-data-files.R`: This script is used sync any data files for individual libraries needed to generate figures to a local folder.
In particular, the `.rds` files for `SCPCS000001`, `SCPCS000216`, and `SCPCS000251` will be stored to folders within the root directory of this repo named `s3_files/SCPCS000001`, `s3_files/SCPCS000216`, and `s3_files/SCPCS000251`, respectively.
Additionally, the results from benchmarking 3 single-cell and 3 single-nuclei samples with both Alevin-fry and Cell Ranger will be stored to a folder within the root directory of this repo named `s3_files/benchmarking_results`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
To run the script use the following command:

```sh
op run -- Rscript sync-data-files.R
```

* `sync-reference-files.R`: This script is used sync any reference files needed to generate figures to a local folder.
These reference files will be stored to a folder within the root directory of this repo named `s3_files/reference_files`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
To run the script use the following command:

```sh
op run -- Rscript sync-reference-files.R
```

* `sync-consensus-celltype-results.R`: This script is used to sync the results from the `cell-type-consensus` module of [`OpenScPCA-nf`](https://github.com/AlexsLemonade/OpenScPCA-nf).
All output files will be saved to a folder within the root directory of this repo named `s3_files/cell-type-consensus-results`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
To run the script you will first need to login to your AWS account with access to `s3://openscpca-nf-workflow-results` and then run the script with the following command:

```sh
Rscript sync-consensus-celltype-results.R --profile <name of AWS profile>
```

* `assign-brain-classifications.R`: This script is used to generate the table assigning each of the brain diagnoses to a subdiagnoses, either High-grade glioma or Low-grade glioma.
The output of this script is `sample-info/brain-classifications-no-multiplexed.tsv` and is used as input to create Figure 5B and Figure 5C.

* `prepare-sample-whitelist.R`: This script is used to prepare the file `../sample-info/sample-whitelist.txt`, which contains a list of samples to include in figures and tables.
This text file is used as input for `Fig1A_sample-disease-barchart.R`, `Fig1B_modality-barchart.R`, and `TableS1_modality-summary.R`.
This script assumes that metadata has been synced with `sync-metadata.R`.
To run the script use the following command:

```sh
Rscript prepare-sample-whitelist.R
```