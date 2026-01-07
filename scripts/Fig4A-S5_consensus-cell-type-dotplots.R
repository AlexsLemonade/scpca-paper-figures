#!/usr/bin/env Rscript

# This script is used to generate the dot plots looking at expression of cell type markers
# in consensus cell types 
# one png file will be saved for each of the 4 diagnoses groups 

library(ggplot2)
library(patchwork)

options(readr.show_col_types = FALSE)

celltype_plotting_functions <- here::here("scripts", "utils", "consensus-celltype-plotting-functions.R")
source(celltype_plotting_functions) # imports `marker_gene_dotplot()`
set.seed(2026) # random marker gene selection; max 10 per validation group

# Set up paths -----------------------------------------------------------------

# output files 
pdf_dir <- here::here("figures", "pdfs")
# Define paths to individual files 
output_pdf_files <- c(
  "Brain and CNS" = file.path(pdf_dir, "Fig4A_brain-dotplot.pdf"),
  "Leukemia" = file.path(pdf_dir, "FigS5A_leukemia-dotplot.pdf"),
  "Sarcoma" = file.path(pdf_dir, "FigS5B_sarcoma-dotplot.pdf"),
  "Other solid tumors" = file.path(pdf_dir, "FigS5C_other-solid-tumors-dotplot.pdf")
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
marker_gene_table_url <- "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/validation-markers.tsv"

# define color palette
color_palette_file <- here::here("palettes", "validation-group-palette.tsv")
celltype_colors <- readr::read_tsv(color_palette_file) |> 
  tibble::deframe()

# Read in metadata -------------------------------------------------------------

# read in project metadata file and filter to project whitelist 
project_whitelist <- readLines(project_whitelist_file)

# read in validation markers as data.tables
markers_df <- readr::read_tsv(marker_gene_table_url) |> 
  # only keep genes unique to a single cell type except HPC which doesn't have any unique genes
  # for HPC we keep all 6 marker genes
  dplyr::filter(gene_observed_count == 1 | validation_group_annotation == "hematopoietic precursor cell") |>
  # keep a random set of 10 marker genes for each validation group
  dplyr::group_by(validation_group_annotation) |>
  dplyr::slice_sample(n = 10) |>
  dplyr::ungroup()


# update celltype_colors: remove colors we don't have marker genes for (but keep "unknown")
keep_colors <- intersect(
  names(celltype_colors), 
  c(unique(markers_df$validation_group_annotation), "unknown")
)
celltype_colors <- celltype_colors[keep_colors]

validation_groups_df <- readr::read_tsv(validation_group_url) |> 
  # rename final assigned group to avoid conflicts when merging in marker gene expression 
  # we want to separate the marker gene group from the actual cell type annotation
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

# Create dot plots -------------------------------------------------------------

# Define width of output PDF for each dotplot
file_widths <- c(24, 16, 24, 24)
file_heights <- c(10, 7, 9, 9)
names(file_widths) <- names(output_pdf_files)
names(file_heights) <- names(output_pdf_files)
dotplot_size_range <- c(5, 6, 5, 5)
names(dotplot_size_range) <- names(output_pdf_files)

plot_list <- output_pdf_files |> 
  purrr::iwalk(\(file, group){
    
    # get only samples in that diag group
    ids <- sample_df |> 
      dplyr::filter(diagnosis_group == group, scpca_sample_id %in% non_multiplex_samples) |> 
      dplyr::pull(scpca_sample_id)
    
    message(glue::glue("Creating dot plot for {group}"))
    
    # create plot
    combined_plot <- marker_gene_dotplot(
      sample_ids = ids,
      consensus_results_dir,
      validation_groups_df,
      markers_df,
      celltype_colors, 
      dotplot_size_range = c(1, dotplot_size_range[group])
    )
    
    # save plot 
    ggsave(file, plot = combined_plot, width = file_widths[group], height = file_heights[group])
    gc() # clean up after each run
  })


