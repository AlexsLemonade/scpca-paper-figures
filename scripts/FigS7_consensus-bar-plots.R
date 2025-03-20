#!/usr/bin/env Rscript

# This script is used to generate the stacked bar chart that shows all cell assignments 
# for all samples in each of the diagnosis groups, except for brain 

library(ggplot2)
library(patchwork)
options(readr.show_col_types = FALSE)

celltype_plotting_functions <- here::here("scripts", "utils", "consensus-celltype-plotting-functions.R")
source(celltype_plotting_functions) # imports `create_celltype_summary()` and `stacked_barchart()`

# Set up paths -----------------------------------------------------------------

# output files 
pdf_dir <- here::here("figures", "pdfs")
output_pdf_files <- c(
  "Leukemia" = file.path(pdf_dir, "FigS7A_leukemia-barplot.pdf"),
  "Sarcoma" = file.path(pdf_dir, "FigS7B_sarcoma-barplot.pdf"),
  "Other solid tumors" = file.path(pdf_dir, "FigS7C_other-solid-tumors-barplot.pdf")
)
# all metadata files 
sample_info_dir <- here::here("sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
diagnosis_groupings_file <- file.path(sample_info_dir, "diagnosis-groupings.tsv")

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

# read in project metadata file and filter to project whitelist 
project_whitelist <- readLines(project_whitelist_file)

validation_groups_df <- readr::read_tsv(validation_group_url) |> 
  # rename final assigned group to match how we do the dot plots and avoid confusion 
  dplyr::select(consensus_annotation, broad_celltype_group = validation_group_annotation)

# read in diagnosis groupings
diagnosis_group_df <- readr::read_tsv(diagnosis_groupings_file) |> 
  dplyr::select(diagnosis = submitted_diagnosis,
                diagnosis_group)

# get sample information
sample_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist) |> 
  dplyr::left_join(diagnosis_group_df, by = c("diagnosis"))

# pull out those that are non-multiplex single cell/nuc
non_multiplex_samples <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(!stringr::str_detect(scpca_sample_id, ";"),
                seq_unit %in% c("cell", "nucleus")) |> 
  dplyr::pull(scpca_sample_id)

# Create bar plots -------------------------------------------------------------

# Define width of output PDF for each barplot
file_widths <- c(33, 20, 20)
names(file_widths) <- names(output_pdf_files)

plot_list <- output_pdf_files |> 
  purrr::iwalk(\(file, group){
    
    # get only samples in that diag group
    ids <- sample_df |> 
      dplyr::filter(diagnosis_group == group, scpca_sample_id %in% non_multiplex_samples) |> 
      dplyr::pull(scpca_sample_id)
    
    message(glue::glue("Creating bar plot for {group}"))
    
    # get list of cell type files using sample ids 
    consensus_results_files <- list.files(
      consensus_results_dir,
      pattern = "_processed_consensus-cell-types\\.tsv\\.gz$",
      recursive = TRUE,
      full.names = TRUE
    )
    celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% ids]
    
    # prep data frame for plotting
    summary_df <- create_celltype_summary(celltype_files, validation_groups_df) |> 
      # join with sample metadata
      dplyr::left_join(sample_df, by = c("project_id" = "scpca_project_id", "sample_id" = "scpca_sample_id"))
    
    # make stacked bar chart
    barplot <- stacked_barchart(summary_df, fill_column = "broad_celltype_group", celltype_colors = celltype_colors)
    
    # save plot 
    ggsave(file, plot = barplot, width = file_widths[[group]], height = 10)
    gc() # clean up after each run
    
  })
