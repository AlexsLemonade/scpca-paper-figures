#!/usr/bin/env Rscript

# This script is used to generate the stacked bar chart that shows all cells in all brain/CNS samples
# samples are faceted based on if they are LGG or HGG 

library(ggplot2)
library(patchwork)
options(readr.show_col_types = FALSE)

celltype_plotting_functions <- here::here("scripts", "utils", "consensus-celltype-plotting-functions.R")
source(celltype_plotting_functions) # imports `stacked_barchart()`

# Set up paths -----------------------------------------------------------------

# output files 
pdf_dir <- here::here("figures", "pdfs")
output_pdf_file <- file.path(pdf_dir, "Fig5B_brain-barchart-all-celltypes.pdf")

# all metadata files 
sample_info_dir <- here::here("sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
brain_classification_file <- file.path(sample_info_dir, "brain-classifications.tsv")

# s3 files 
s3_dir <- here::here("s3_files")
sample_metadata_file <- file.path(s3_dir, "scpca-sample-metadata.tsv")
library_metadata_file <- file.path(s3_dir, "scpca-library-metadata.tsv")
consensus_results_dir <- file.path(s3_dir, "cell-type-consensus-results")

# validation groups and marker gene table urls 
validation_group_url <- "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/consensus-validation-groups.tsv"

# define color palette
color_palette_file <- here::here("palettes", "validation-group-palette.tsv")
celltype_colors <- readr::read_tsv(color_palette_file) |> 
  tibble::deframe()

# Prep metadata ----------------------------------------------------------------

# get validation groups to use for grouping cells 
validation_groups_df <- readr::read_tsv(validation_group_url) |> 
  # rename final assigned group to match how we do the dot plots and avoid confusion 
  dplyr::select(consensus_annotation, broad_celltype_group = validation_group_annotation)

# brain projects 
brain_project_ids <- c("SCPCP000001", "SCPCP000002", "SCPCP000010", "SCPCP000021", "SCPCP000009")
# pull out those that are non-multiplex single cell/nuc
non_multiplex_samples <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(!stringr::str_detect(scpca_sample_id, ";"),
                seq_unit %in% c("cell", "nucleus")) |> 
  dplyr::pull(scpca_sample_id)

# brain classifications
brain_classification_df <- readr::read_tsv(brain_classification_file) |> 
  dplyr::select(diagnosis = submitted_diagnosis, subdiagnosis_group)

# read in sample metadata and select samples
sample_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% brain_project_ids & scpca_sample_id %in% non_multiplex_samples) |> 
  # add in classifications
  dplyr::left_join(brain_classification_df, by = "diagnosis")

sample_ids <- sample_df$scpca_sample_id

# Prep and plot cell types -----------------------------------------------------

# get list of cell type files using sample ids 
consensus_results_files <- list.files(
  consensus_results_dir,
  pattern = "_processed_consensus-cell-types\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)
celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]

# prep data frame for plotting
summary_df <- create_celltype_summary(celltype_files, validation_groups_df) |> 
  # join with sample metadata
  dplyr::left_join(sample_df, by = c("project_id" = "scpca_project_id", "sample_id" = "scpca_sample_id")) |> 
  # only keep HGG and LGG 
  dplyr::filter(subdiagnosis_group %in% c("High-grade glioma", "Low-grade glioma")) 

# make stacked bar chart and facet 
hgg_lgg_barplot <- stacked_barchart(summary_df, fill_column = "broad_celltype_group", celltype_colors = celltype_colors, facet_variable = "subdiagnosis_group")

# save plot 
ggsave(output_pdf_file, plot = hgg_lgg_barplot, width = 20, height = 10)

