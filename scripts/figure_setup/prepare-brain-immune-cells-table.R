#!/usr/bin/env Rscript

# This script is used to create a table with immune cell types 
# to be displayed in the dotplot for Fig 4D
# we also include the validation group for each immune cell type as those are used to get marker genes to display 

options(readr.show_col_types = FALSE)

# Set up -----------------------------------------------------------------------

# define the cell type order that we want in the dotplot
celltype_order <- c(
  "B cell", 
  "lymphocyte of B lineage",
  "naive B cell",
  "mature B cell",
  "memory B cell",
  "plasma cell",
  "T cell",
  "mature T cell",
  "mature alpha-beta T cell",
  "CD8-positive, alpha-beta T cell",
  "CD8-positive, alpha-beta memory T cell",
  "CD4-positive, alpha-beta T cell",
  "memory T cell",
  "regulatory T cell",
  "mononuclear phagocyte",
  "myeloid leukocyte",
  "macrophage",
  "tissue-resident macrophage",
  "monocyte",
  "microglial cell",
  "dendritic cell",
  "natural killer cell"
)

# all metadata files 
sample_info_dir <- here::here("sample-info")

# s3 files 
s3_dir <- here::here("s3_files")
sample_metadata_file <- file.path(s3_dir, "scpca-sample-metadata.tsv")
library_metadata_file <- file.path(s3_dir, "scpca-library-metadata.tsv")
consensus_results_dir <- file.path(s3_dir, "cell-type-consensus-results")

# validation group url 
validation_group_url <- "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/consensus-validation-groups.tsv"

# output file 
immune_cells_file <- file.path(sample_info_dir, "brain-immune-celltypes.tsv")

# Get consensus results for brain samples --------------------------------------

# brain projects 
brain_project_ids <- c("SCPCP000001", "SCPCP000002", "SCPCP000010", "SCPCP000021", "SCPCP000009")
# pull out those that are non-multiplex single cell/nuc
non_multiplex_samples <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(!stringr::str_detect(scpca_sample_id, ";"),
                seq_unit %in% c("cell", "nucleus")) |> 
  dplyr::pull(scpca_sample_id)

# get sample information
# read in sample metadata and select samples to include
sample_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% brain_project_ids & scpca_sample_id %in% non_multiplex_samples)
sample_ids <- sample_df$scpca_sample_id

# get a list of the consensus results files for specified samples
consensus_results_files <- list.files(
  consensus_results_dir,
  pattern = "_consensus-cell-types\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)
celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]

# read in all consensus results 
consensus_df <- duckplyr::read_csv_duckdb(celltype_files, options = list(sep = "\t", union_by_name = TRUE))

# Check for any missing consensus cell types -----------------------------------

# all consensus cell types present with > 50 cells 
brain_consensus_celltypes <- consensus_df |> 
  dplyr::count(consensus_annotation) |> 
  dplyr::filter(n > 50) |> 
  dplyr::pull(consensus_annotation)

# list of any that will not be included in the dot plot
non_immune_celltypes <- setdiff(brain_consensus_celltypes, celltype_order) |> 
  paste0(collapse = "\n")

# print them out so we can check that this is correct
message(
  glue::glue("The following cell types are found but will not be showin the dotplot:
             
             {non_immune_celltypes}
             ")
)

# Prepare immune cell reference table ------------------------------------------

# read in existing validation groups
validation_group_df <- readr::read_tsv(validation_group_url)

# export immune cell order as a data frame
immune_cells_df <- data.frame(
  consensus_annotation = celltype_order
) |> 
  dplyr::left_join(validation_group_df, by = c("consensus_annotation"))

readr::write_tsv(immune_cells_df, "sample-info/brain-immune-celltypes.tsv")

