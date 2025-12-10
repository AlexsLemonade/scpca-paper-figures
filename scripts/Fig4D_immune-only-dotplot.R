#!/usr/bin/env Rscript

# This script is used to generate the dot plots looking at expression of cell type markers
# in immune consensus cell types for brain and CNS tumors included in Fig 4

library(ggplot2)
library(patchwork)

options(readr.show_col_types = FALSE)

celltype_plotting_functions <- here::here("scripts", "utils", "consensus-celltype-plotting-functions.R")
source(celltype_plotting_functions) # imports `make_dotplot()`

# Define file paths ------------------------------------------------------------

# all metadata files 
sample_info_dir <- here::here("sample-info")

# s3 files 
s3_dir <- here::here("s3_files")
sample_metadata_file <- file.path(s3_dir, "scpca-sample-metadata.tsv")
library_metadata_file <- file.path(s3_dir, "scpca-library-metadata.tsv")
consensus_results_dir <- file.path(s3_dir, "cell-type-consensus-results")

# color palette 
palette_file <- here::here("palettes", "immune-palette.tsv")
immune_colors <- readr::read_tsv(palette_file) |> 
  tibble::deframe()

# marker gene table url
marker_gene_table_url <- "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/consensus-markers.tsv"

# all possible consensus immune cell types in the brain and their associated lineages 
immune_cells_file <- file.path(sample_info_dir, "brain-immune-celltypes.tsv")

# output file 
output_pdf_file <- here::here("figures", "pdfs", "Fig4D_brain-immune-dotplot.pdf")

# Prep cell type and marker gene info ------------------------------------------


# table of immune cell types and lineages
immune_cells_df <- readr::read_tsv(immune_cells_file) |> 
  dplyr::select(consensus_annotation, marker_gene_lineage)

# define the order of celltypes to use throughout 
celltype_order <- readr::read_tsv(immune_cells_file) |>
  dplyr::pull(consensus_annotation)

# read in validation markers
markers_df <- readr::read_tsv(marker_gene_table_url) |> 
  #dplyr::filter(consensus_annotation %in% all_immune_cells) |> 
  dplyr::right_join(immune_cells_df, by = c("consensus_annotation")) |> 
  tidyr::drop_na() |> # drop cell types with no marker genes 
  dplyr::rename("validation_group_annotation" = consensus_annotation)

  # only keep markers that are in our immune cell list
  dplyr::filter(consensus_annotation %in% celltype_order) |> 
  dplyr::rename("validation_group_annotation" = consensus_annotation) 

# get marker genes that are unique to the included cell types
unique_markers <- markers_df |> 
  dplyr::count(gene_symbol) |> 
  # only keep genes that are unique to the cell types being shown
  dplyr::filter(n <= 1) |> 
  dplyr::pull(gene_symbol)

lineage_markers <- markers_df |> 
  dplyr::select(gene_symbol, marker_gene_lineage) |> 
  unique() |> 
  dplyr::count(gene_symbol) |> 
  dplyr::filter( n <= 1) |> 
  dplyr::pull(gene_symbol)

# filter to genes that are only found in one lineage 
markers_df <- markers_df |> 
  dplyr::filter(gene_symbol %in% unique_markers)

# Get brain sample ids ---------------------------------------------------------

# brain projects 
brain_project_ids <- c("SCPCP000001", "SCPCP000002", "SCPCP000010", "SCPCP000021", "SCPCP000009")
# pull out those that are non-multiplex single cell/nuc
non_multiplex_samples <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(!stringr::str_detect(scpca_sample_id, ";"),
                seq_unit %in% c("cell", "nucleus")) |> 
  dplyr::pull(scpca_sample_id)


# get sample information
# read in sample metadata and select samples
sample_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% brain_project_ids & scpca_sample_id %in% non_multiplex_samples)


sample_ids <- sample_df$scpca_sample_id

# Read in consensus and gene exp files -----------------------------------------

# This code is modified from the `marker_gene_dotplot()` function
# list of consensus files 
consensus_results_files <- list.files(
  consensus_results_dir,
  pattern = "_consensus-cell-types\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)
celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]

# gene expression used for dot plots
gene_exp_files <- list.files(
  consensus_results_dir,
  pattern = "_marker-gene-expression\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)
gene_exp_files <- gene_exp_files[basename(dirname(gene_exp_files)) %in% sample_ids]

# Use DuckDB for dplyr functions
suppressMessages(duckplyr::methods_overwrite())

# convert to duckdb (intermediate data frame to avoid duckdb error with readr objects)
markers_df <- as.data.frame(markers_df) |> duckplyr::as_duckdb_tibble()

consensus_df <- duckplyr::read_csv_duckdb(celltype_files, options = list(sep = "\t", union_by_name = TRUE)) |> 
  dplyr::filter(consensus_annotation %in% celltype_order)

gene_exp_df <- duckplyr::read_csv_duckdb(gene_exp_files, options = list(sep = "\t", union_by_name = TRUE)) |>
  dplyr::mutate(detected = logcounts > 0) 

# Join all consensus results and marker gene info
consensus_df <- consensus_df |>
  dplyr::left_join(gene_exp_df, by = c("barcodes", "library_id")) |>
  # add marker gene information (associated validation group annotation, gene observed count, percent tissues)
  # account for the same gene being present in multiple cell types
  dplyr::left_join(markers_df, by = c("ensembl_gene_id"), relationship = "many-to-many")

# Prep for dotplot -------------------------------------------------------------
# get total number of cells per final annotation group
total_cells_df <- consensus_df |>
  dplyr::select(library_id, barcodes, consensus_annotation) |>
  dplyr::distinct() |>
  dplyr::count(consensus_annotation, name = "total_cells")

# table with one row per unique broad cell type/ marker gene combination
# first all cells in with the same consensus_annotation are grouped together
# then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
# do this for every possible marker gene/ validation group annotation assignment

group_stats_df <- consensus_df |>
  # for each assigned cell type/marker gene combo get total detected and mean expression
  # group by both consensus and validation group to account for genes that are expressed in more than one cell type
  # here the "validation group" is just the consensus cell type associated with the marker gene
  # the consensus annotation column is the assignment associated with that cell
  dplyr::summarize(
    .by = c("consensus_annotation", "ensembl_gene_id", "validation_group_annotation"),
    detected_count = sum(detected),
    mean_exp = mean(logcounts)
  ) |>
  # add in validation group for marker genes
  # this includes all possible marker genes and all possible validation group assignments
  dplyr::left_join(markers_df, by = c("ensembl_gene_id", "validation_group_annotation"), relationship = "many-to-many") |> # add total cells
  dplyr::left_join(total_cells_df, by = c("consensus_annotation")) |>
  # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild
  dplyr::filter(total_cells > 50) |>
  dplyr::mutate(
    # get total percent expressed
    percent_exp = (detected_count / total_cells) * 100,
    # account for NA/unknowns and set axes order
    consensus_annotation = factor(consensus_annotation, levels = c(celltype_order, "unknown"))
  )

# get list of celltypes to keep for plotting
celltype_groups <- group_stats_df |>
  dplyr::pull(consensus_annotation) |>
  unique() |>
  as.character()

# filter markers to those that are actually relevant
# we will only plot the marker genes for cell types that are present and only include marker genes that are present
filtered_markers_df <- markers_df |>
  dplyr::filter(
    validation_group_annotation %in% celltype_groups & gene_symbol %in% group_stats_df$gene_symbol
  ) |> 
  # ensure order matches the order of the legend
  dplyr::mutate(
    validation_group_annotation = factor(validation_group_annotation, levels = celltype_order),
    # account for marker genes being present in more than one cell type
    gene_symbol = make.unique(gene_symbol),
  ) |> 
  dplyr::arrange(validation_group_annotation)

# specify x axis order for dotplot
marker_gene_order <- filtered_markers_df |>
  dplyr::pull(gene_symbol)

# filter out low expressed genes
dotplot_df <- group_stats_df |>
  dplyr::filter(mean_exp > 0, percent_exp > 10) |>
  dplyr::arrange(consensus_annotation) |>
  # add a label for the plot
  dplyr::mutate(y_label = glue::glue("{consensus_annotation} ({total_cells})")) |>
  # remove marker genes that aren't present in final annotations and set x axis order
  dplyr::filter(gene_symbol %in% marker_gene_order) |>
  dplyr::mutate(
    # set orders of gene symbol and validation groups
    y_label = factor(y_label, levels = rev(unique(y_label))),
    gene_symbol = factor(gene_symbol, levels = marker_gene_order),
    validation_group_annotation = factor(validation_group_annotation, levels = celltype_order)
  )

# Make dotplot -----------------------------------------------------------------

combined_plot <- make_dotplot(
  dotplot_df,
  cell_type_label = "Consensus cell type annotation", 
  celltype_colors = immune_colors
)

ggsave(output_pdf_file, plot = combined_plot, width = 30, height = 10)

suppressMessages(duckplyr::methods_restore()) # back to dplyr
