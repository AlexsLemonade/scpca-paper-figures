# Scripts for generating figures and tables

This folder contains all scripts used to generate figures and tables.
Before running these scripts, you will need to prepare data for script input:

* If you are a Data Lab staff member, you will need to run scripts described in the **Figure set up** section at the bottom of this document
* If you are not a Data Lab staff member, you will need to follow instructions provided in the file `../reproduce-figures/obtain-prepare-data.md`

## Generating figures and tables

The following scripts can be used to generate figures and tables:

* `Fig1A_sample-disease-barchart.R`: This script is used to generate Figure 1A, which includes a summary of the types of diagnoses found in the ScPCA portal.

* `Fig1B_modality-barchart.R`: This script is used to generate Figure 1B, which includes a summary of the types of modalities found in the ScPCA portal.

* `Fig2B-F_qc-plots.R`: This script is used to generate panels Figure 2B-F, which includes simplified and miniature versions of the plots found in the main QC report included with each sample download.

* `Fig3C_singler-cellassign-heatmap.R`: This script is used to generate figure 3C, which is a heatmap comparing cell type annotations from automated methods to consensus cell types.

* `Fig4A-S5_consensus-cell-type-dotplots.R`: This script is used to generate Figure 4A and supplemental figure 5A-C.
Figure 4A includes a dot plot summarizing marker gene expression across all consensus cell types in the Brain and CNS samples.
Supplemental Figure 5A-C include dot plots summarizing marker gene expression across all consensus cell types in the other diagnosis groups, leukemia, sarcoma, and other solid tumors.

* `Fig4B-C_brain-barplots.R`: This script is used to generate Figures 4B and 4C.
Figure 4B shows the percentage of each library that is annotated as each of the consensus cell types for all libraries in High-grade and Low-grade glioma samples.
Figure 4C shows the percentage of each library that is annotated as each of the _immune_ consensus cell types, emphasizing myeloid and lymphocyte cell types, for all libraries in High-grade and Low-grade glioma samples.

* `Fig4D_immune-only-dotplot.R`: This script is used to generate Figure 4D, which shows the expression of marker genes for all immune consensus cell types in Brain and CNS samples in a dot plot.

* `Fig5A-B_nb-umap-panels.R`: This script is used to generate Figures 5A and 5B.
Figure 5A shows an unintegrated UMAP of Neuroblastoma samples from project `SCPCP000004` colored by broad cell type annotation as performed through the OpenScPCA Project.
Figure 5B shows an unintegrated UMAP of Neuroblastoma samples from project `SCPCP000004` colored by the total CNV as inferred by `inferCNV`.

* `Fig5C_nb-cnv-heatmaps.R`: This script is used to generate Figure 5C, which contains CNV heatmaps for chromosomes with canonical events for a single Neuroblastoma sample.

* `Fig5D_nb-ridgeplot.R`: This script is used to generate Figure 5D, which contains a ridgeplot of CNV events across top consensus cell types for a single Neuroblastoma sample.

* `Fig6-FigS7_bulk-analysis.R `: This script is used to generate Figures 6 and S8.
Figure 6 shows (A) scatterplots of the relationship between bulk and pseudobulk counts for relevant Brain and CNS projects, and (B) bar plots of odds ratios from overrepresentation analysis of bulk expression data for relevant Brain and CNS projects.
Figure S7 shows (A) scatterplots of the relationship between bulk and pseudobulk counts for projects not shown in Figure 6A, and (B) bar plots of odds ratios from overrepresentation analysis of bulk expression data for projects not shown in Figure 6B.

* `FigS1A_memory-time-comparison.R`: This script is used to generate supplemental Figure 1A, which shows a comparison of total run time and peak memory usage for Cell Ranger and Alevin-fry.
This script uses the trace files found in `nextflow_logs/`.

* `FigS1B-D_method-metrics-comparison.R`: This script is used to generate supplemental figures 1B-D, which compares cell and gene level metrics between libraries quantified using Cell Ranger and Alevin-fry.

* `FigS2B-D_adt-plots.R`: This script is used to generate supplemental Figure 2B-D, which includes simplified and miniature versions of the plots in the ADT section of the main QC report included with each sample download.

* `FigS3D_merged-umaps.R`: This script is used to generate Figure 3SD, which includes an example of the UMAPs shown in the merged report.

* `FigS4_t-cell-umap.R`: This script is used to generate Supplemental figure 4A, which shows a faceted UMAP of the top T cell types labeled by each automated method for `SCPCL000049`. 

* [TODO] `FigS4_submitter-celltypes-heatmap.R`: This script is used to generate supplemental Figure 5, which includes a heatmap comparing the `CellAssign` and `SingleR` annotations to submitter-provided annotations for an example library.

* `FigS6_consensus-bar-plots.R`: This script is used to generate Supplemental figure 7A-C, which show the percentage of each library that is annotated as each of the consensus cell types for all libraries in each diagnosis group, excluding the Brain and CNS samples.

* `TableS1_modality-summary.R`: This script is used to generate supplemental Table 1, which contains a summary of the types of libraries found in each project.

* `TableS2_cellassign-ref-summary.R`: This script is used to generate supplemental Table 2, which includes one row for each ScPCA project on the Portal and the associated diagnoses, reference used from `PanglaoDB`, and list of organs used to construct the reference.

## Old figures

The `old_figures` folder contains the scripts used to generate figures present in previous versions of the manuscript.

* `FigS4_celldex-ref-comparison.R `: This script is used to generate a previous version of supplemental Figure 4, which included a figure comparing `celldex` references across samples.

* `FigS5A_cellassign-justification.R`: This script is used to generate a previous version of supplemental Figure 5A, which included a faceted UMAP of automated cell type annotations determined by `CellAssign`.

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

* `prepare-brain-immune-cells-table.R`: This script is used to prepare the file `../sample-info/brain-immune-celltypes.tsv`, which contains a table of all immune consensus cell types with at least 50 assigned cells in the brain and CNS samples. 
This table is used to specify which cell types should be shown in the dot plot displaying all immune cell types and their associated marker genes.
