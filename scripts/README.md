# Scripts for generating figures and tables

This folder contains all scripts used to generate figures and tables.

## Figure set up

The `figure_setup` folder contains scripts that are required to be run _prior to_ generating figures.

1. `sync-metadata.R`: This script is used to sync any metadata files found on S3 to a local folder.
In particular, the `scpca-sample-metadata.tsv` and `scpca-library-metadata.tsv` will be stored to a folder within the root directory of this repo named `s3_files`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
The `s3_files` folder is included in the `.gitignore`, so the first time you go to generate figures you will need to run this script as these files are not available in the repo.
To run the script use the following command:

```sh
op run -- Rscript sync-metadata.R
```

2. `sync-data-files.R`: This script is used sync any data files for individual libraries needed to generate figures to a local folder.
In particular, the `.rds` files for `SCPCS000001`, `SCPCS000216`, and `SCPCS000251` will be stored to folders within the root directory of this repo named `s3_files/SCPCS000001`, `s3_files/SCPCS000216`, and `s3_files/SCPCS000251`, respectively.
Additionally, the results from benchmarking 3 single-cell and 3 single-nuclei samples with both Alevin-fry and Cell Ranger will be stored to a folder within the root directory of this repo named `s3_files/benchmarking_results`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
To run the script use the following command:

```sh
op run -- Rscript sync-data-files.R
```

3. `sync-reference-files.R`: This script is used sync any reference files needed to generate figures to a local folder.
These reference files will be stored to a folder within the root directory of this repo named `s3_files/reference_files`.
In order to generate some of the figures (see more on which figures require this script below), this script will be need to run first.
To run the script use the following command:

```sh
op run -- Rscript sync-reference-files.R
```

## Generating figures and tables

The following scripts can be used to generate figures and tables:

1. `Fig1A_sample-disease-barchart.R`: This script is used to generate Figure 1A, which includes a summary of the types of diagnoses found in the ScPCA portal.
Before running this script, you must run `figure_setup/sync-metadata.R`.

2. `Fig1B_modality-barchart.R`: This script is used to generate Figure 1B, which includes a summary of the types of modalities found in the ScPCA portal.
Before running this script, you must run `figure_setup/sync-metadata.R`.

3. `TableS1_modality-summary.R`: This script is used to generate supplemental Table 1, which contains a summary of the types of libraries found in each project.
Before running this script, you must run `figure_setup/sync-metadata.R`.

4. `Fig2B-G_qc-plots.R`: This script is used to generate panels Figure 2B-G, which includes simplified and miniature versions of the plots found in the main QC report included with each sample download.
Before running this script, you must run `figure_setup/sync-data-files.R`.

5. `FigS1A_memory-time-comparison.R`: This script is used to generate supplemental Figure 1A, which shows a comparison of total run time and peak memory usage for Cell Ranger and Alevin-fry.
This script uses the trace files found in `nextflow_logs`.

6. `FigS1B-D_method-metrics-comparison.R`: This script is used to generate supplemental figures 1B-D, which compares cell and gene level metrics between libraries quantified using Cell Ranger and Alevin-fry.
Before running this script, you must run `figure_setup/sync-data-files.R` and `sync-reference-files.R`.

7. `FigS2B-D_adt-plots.R`: This script is used to generate supplemental Figure 2B-D, which includes simplified and miniature versions of the plots in the ADT section of the main QC report included with each sample download.
Before running this script, you must run `figure_setup/sync-data-files.R`.

8. `Fig3D_merged-umaps.R`: This script is used to generate Figure 3D, which includes an example of the UMAPs shown in the merged report.
Before running this script, you must run `figure_setup/sync-data-files.R`.

9. `Fig4B_singler-cellassign-heatmap.R`: This script is used to generate Figure 4B, which includes an example of the heatmap comparing cell type annotations from `SingleR` to `CellAssign`.
Before running this script, you must run `figure_setup/sync-data-files.R`.

10. `TableS2_cellassign-ref-summary.R`: This script is used to generate supplemental Table 2, which includes one row for each ScPCA project on the Portal and the associated diagnoses, reference used from `PanglaoDB`, and list of organs used to construct the reference.
Before running this script, you must run `figure_setup/sync-metadata.R`.

11. `FigS4_celldex-ref-comparison.R`: This script is used to generate supplemental Figure 6, which compares the delta median statistic calculated from running `SingleR` on a subset of ScPCA libraries with 4 different `celldex` references.
Before running this script, you must run both `figure_setup/sync-data-files.R` and `figure_setup/sync-metadata.R`.

12. `FigS5A-B_cellassign-justification.R`: This script is used to generate supplemental Figure 7A-B, which includes a UMAP summarizing cell type annotations using `CellAssign` and a heatmap comparing the `CellAssign` annotations to submitter-provided annotations.
Before running this script, you must run `figure_setup/sync-data-files.R`.

13. `FigS6A-B_celltype-diagnostic-plots.R`: This script is used to generate supplemental Figure 4A-B, which includes examples of the diagnostic plots for annotating cell types with `SingleR` and `CellAssign`.
Before running this script, you must run `figure_setup/sync-data-files.R`.

14. `FigS7_submitter-celltypes-heatmap.R`: This script is used to generate supplemental Figure 5, which includes an example heatmap comparing submitter provided annotations to automated annotations from `SingleR` and `CellAssign`.
Before running this script, you must run `figure_setup/sync-data-files.R`.
