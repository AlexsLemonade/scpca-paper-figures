#!/usr/bin/env Rscript

# This file contains helper functions for plotting consensus cell types

# packages required for functions
library(ggplot2)
library(patchwork)
library(data.table)
library(duckplyr)

#' Dot plot showing expression of marker genes across assigned cell types
#'
#' @param sample_ids Character vector of ScPCA sample ids to include in plot
#' @param consensus_results_dir Directory where results from cell-type-consensus module of OpenScPCA-nf lives
#' @param validation_groups_df Data frame assigning consensus cell types to broader validation groups
#' @param markers_df Data frame with marker genes for each cell type
#' @param celltype_colors Named vector of colors to use for each broader validation group
#'
#' @returns Dot plot with summarized expression of marker genes for consensus cell types
marker_gene_dotplot <- function(
    sample_ids,
    consensus_results_dir,
    validation_groups_df,
    markers_df,
    celltype_colors) {
  # list all cell type assignments files
  consensus_results_files <- list.files(
    consensus_results_dir,
    pattern = "_processed_consensus-cell-types\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]

  # gene expression used for dot plots
  gene_exp_files <- list.files(
    consensus_results_dir,
    pattern = "_processed_marker-gene-expression\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  gene_exp_files <- gene_exp_files[basename(dirname(gene_exp_files)) %in% sample_ids]

  # convert to duckdb (intermediate tibble to avoid duckdb error)
  validation_groups_df <- as_tibble(validation_groups_df) |> as_duckdb_tibble()
  markers_df <- as_tibble(markers_df) |> as_duckdb_tibble()

  # read in files
  consensus_df <- celltype_files |>
    purrr::map(fread) |>
    data.table::rbindlist(fill = TRUE, use.names = TRUE) |>
    as_duckdb_tibble(prudence = "thrifty")


  # use a tempfile and read_csv_duckdb() to reduce memory usage for reading in gene expression files
  gene_exp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(gene_exp_csv), add = TRUE)
  csv_append <- FALSE
  for (f in gene_exp_files) {
    gene_exp_df <- fread(f)
    fwrite(gene_exp_df, gene_exp_csv, append = csv_append)
    csv_append <- TRUE # after the first, append
  }
  gene_exp_df <- read_csv_duckdb(gene_exp_csv) |>
    mutate(detected = logcounts > 0)

  # Join all consensus results and marker gene info
  consensus_df <- consensus_df |>
    # add in broad cell type group which is used for plotting
    # groups similar cell types together
    left_join(validation_groups_df, by = "consensus_annotation") |>
    left_join(gene_exp_df, by = c("barcodes", "library_id")) |>
    # add marker gene information (associated validation group annotation, gene observed count, percent tissues)
    # account for the same gene being present in multiple cell types
    left_join(markers_df, by = "ensembl_gene_id", relationship = "many-to-many")

  # prep for plots
  # get total number of cells per final annotation group
  total_cells_df <- consensus_df |>
    select(library_id, barcodes, broad_celltype_group) |>
    distinct() |>
    count(broad_celltype_group, name = "total_cells")

  # table with one row per unique broad cell type/ marker gene combination
  # first all cells in with the same broad_celltype_group (determined based on consensus_annotation) are grouped together
  # then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
  # do this for every possible marker gene/ validation group assignment
  # second we calculate the mean expression and mean percentage of all marker genes in a given validation group (this value is used only in the second section of the report)
  group_stats_df <- consensus_df |>
    # for each assigned cell type/marker gene combo get total detected and mean expression
    # group by both broad group and validation group to account for genes that are expressed in more than one cell type
    summarize(
      .by = c("broad_celltype_group", "ensembl_gene_id", "validation_group_annotation"),
      detected_count = sum(detected),
      mean_exp = mean(logcounts)
    ) |>
    # add in validation group for marker genes
    # this includes all possible marker genes and all possible validation group assignments
    left_join(markers_df, by = c("ensembl_gene_id", "validation_group_annotation"), relationship = "many-to-many") |>
    # now get the mean expression/ mean percentage across all marker genes for a given validation group
    # here the broad_celltype_group is the final assigned annotation for that group of cells
    # the validation_group_annotation refers to the cell type that marker gene is associated with
    mutate(
      .by = c("broad_celltype_group", "validation_group_annotation"),
      # calculate mean expression/detected across all markers for a specific group
      all_markers_mean_exp = mean(mean_exp),
      all_markers_detected_count = mean(detected_count)
    ) |> # add total cells
    left_join(total_cells_df, by = c("broad_celltype_group")) |>
    # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild
    filter(total_cells > 50) |>
    mutate(
      # get total percent expressed
      percent_exp = (detected_count / total_cells) * 100,
      all_markers_percent_exp = (all_markers_detected_count / total_cells) * 100,
      # account for NA/unknowns and set axes order
      broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown") |>
        factor(levels = c(unique(markers_df$validation_group_annotation), "unknown"))
    )

  # get list of celltypes to keep and assign colors
  celltype_groups <- group_stats_df |>
    pull(broad_celltype_group) |>
    unique() |>
    as.character()

  # filter markers to those that are actually relevant
  # we will only plot the marker genes for cell types that are part of the assigned broad validation group for this group of samples
  # we don't care about plotting marker genes for cell types that aren't present here
  filtered_markers_df <- markers_df |>
    filter(
      validation_group_annotation %in% celltype_groups,
      gene_symbol %in% group_stats_df$gene_symbol
    )

  # specify x axis order for dotplot
  marker_gene_order <- filtered_markers_df |>
    pull(gene_symbol)

  # set order for cell types
  celltype_order <- unique(filtered_markers_df$validation_group_annotation)

  # filter out low expressed genes
  dotplot_df <- group_stats_df |>
    filter(mean_exp > 0, percent_exp > 10) |>
    arrange(broad_celltype_group) |>
    # add a label for the plot
    mutate(y_label = as.factor(glue::glue("{broad_celltype_group} ({total_cells})"))) |>
    # remove marker genes that aren't present in final annotations and set x axis order
    filter(gene_symbol %in% marker_gene_order) |>
    mutate(
      # set orders of gene symbol and validation groups
      gene_symbol = factor(gene_symbol, levels = marker_gene_order),
      validation_group_annotation = factor(validation_group_annotation, levels = celltype_order)
    ) |>
    as_tibble()


  # make dotplot with marker gene exp
  dotplot <- ggplot(dotplot_df, aes(y = forcats::fct_rev(y_label), x = gene_symbol, color = mean_exp, size = percent_exp)) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    facet_grid(cols = vars(validation_group_annotation), scales = "free", space = "free") +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "transparent", color = NA),
      strip.placement = "outside",
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.ticks.x = element_blank(),
      text = element_text(size = 14),
      panel.spacing = unit(0.5, "lines") # adjust spacing and match with annotation bar
    ) +
    labs(
      x = "",
      y = "Broad cell type annotation",
      color = "Mean gene expression",
      size = "Percent cells expressed"
    )


  # add annotation bar aligning marker genes with validation group
  color_bar <- ggplot(dotplot_df, aes(x = gene_symbol, y = 1, fill = validation_group_annotation)) +
    geom_tile() +
    facet_grid(cols = vars(validation_group_annotation), scales = "free", space = "free") +
    scale_fill_manual(values = celltype_colors, breaks = levels(dotplot_df$validation_group_annotation)) +
    ggmap::theme_nothing() +
    theme(
      strip.background = element_rect(fill = "transparent", color = NA),
      strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 12),
      strip.placement = "outside",
      legend.position = "none",
      panel.spacing = unit(0.5, "lines"),
      strip.clip = "off"
    ) +
    labs(fill = "")

  combined_plot <- color_bar / dotplot +
    patchwork::plot_layout(ncol = 1, heights = c(0.1, 4))

  return(combined_plot)
}
