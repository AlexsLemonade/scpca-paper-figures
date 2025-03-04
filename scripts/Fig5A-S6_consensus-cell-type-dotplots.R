#!/usr/bin/env Rscript

# This script is used to generate the dot plots looking at expression of cell type markers
# in consensus cell types 
# one png file will be saved for each of the 4 diagnoses groups 

# load project
renv::load()

library(ggplot2)
library(patchwork)
library(data.table)

# Set up paths -----------------------------------------------------------------

celltype_plotting_functions <- here::here("scripts", "utils", "consensus-celltype-plotting-functions.R")
source(celltype_plotting_functions)

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

# color palette
# TODO: Turn this into a TSV that's saved in the repo 
celltype_colors <- c(
  # lymphocytes
  "B cell" = "#AA0DFE",
  "plasma cell" = "#782AB6",
  "T cell" = "#3283FE",
  "innate lymphoid cell" = "#325A9B",
  # myeloid
  "dendritic cell" = "#16FF32",
  "macrophage" = "#1C8356",
  "monocyte" = "#1CBE4F",
  "myeloid" = "#90AD1C",
  "natural killer cell" = "#FBE426",
  # hpsc 
  "hematopoietic precursor cell" = "#FEAF16",
  # stem cell
  "stem cell" = "#F7E1A0",
  # stromal cells
  "adipocyte" = "#FA0087",
  "chondrocyte" = "#B00068",
  "endothelial cell" = "#FC1CBF",
  "epithelial cell" = "#C075A6",
  "fibroblast" = "#B10DA1",
  "melanocyte" = "#FE00FA",
  "muscle cell" = "#F6222E",
  "pericyte" = "#F8A19F",
  "stromal cell" = "#C4451C",
  # neural cells
  "neuron" = "#1CFFCE",
  "astrocyte" = "#2ED9FF",
  "glial cell" = "#565656",
  # other 
  "mesangial cell" = "#E2E2E2"
)

# output files 
pdf_dir <- here::here("figures", "pdfs")
# Define paths to individual files 
output_pdf_files <- c(
  "Brain and CNS" = file.path(pdf_dir, "Fig5A_brain-dotplot.pdf"),
  #"Leukemia" = file.path(pdf_dir, "FigS6A_leukemia-dotplot.pdf"),
  "Sarcoma" = file.path(pdf_dir, "FigS6B_sarcoma-dotplot.pdf"),
  "Other solid tumors" = file.path(pdf_dir, "FigS6C_other-solid-tumors-dotplot.pdf")
)

# Read in metadata -------------------------------------------------------------

# read in project metadata file and filter to project whitelist 
project_whitelist <- readLines(project_whitelist_file)

# read in validation markers as data.tables
markers_dt <- fread(marker_gene_table_url) |> 
  # only keep genes unique to a single cell type except HPC which doesn't have any unique genes
  # for HPC we keep all 6 marker genes
  dplyr::filter(gene_observed_count == 1 | validation_group_annotation == "hematopoietic precursor cell")

validation_groups_dt <- fread(validation_group_url) |> 
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
      validation_groups_dt,
      markers_dt,
      celltype_colors
    )
    
    # save plot 
    ggsave(file, plot = combined_plot, width = 18, height = 10)
    
  })


