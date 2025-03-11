#!/usr/bin/env Rscript

# This file contains helper functions for plotting consensus cell types

# packages required for functions
library(ggplot2)
library(patchwork)

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

  # Use DuckDB for dplyr functions
  suppressMessages(duckplyr::methods_overwrite())
  
  # convert to duckdb (intermediate data frame to avoid duckdb error with readr objects)
  validation_groups_df <- as.data.frame(validation_groups_df) |> duckplyr::as_duckdb_tibble()
  markers_df <- as.data.frame(markers_df) |> duckplyr::as_duckdb_tibble()

  # read in files directly to duckdb tables
  # specify all varchar for consensus to avoid parsing error
  consensus_df <- duckplyr::read_csv_duckdb(celltype_files, options = list(sep = "\t", union_by_name = TRUE)) 

  gene_exp_df <- duckplyr::read_csv_duckdb(gene_exp_files, options = list(sep = "\t", union_by_name = TRUE)) |>
    dplyr::mutate(detected = logcounts > 0)

  # Join all consensus results and marker gene info
  consensus_df <- consensus_df |>
    # add in broad cell type group which is used for plotting
    # groups similar cell types together
    dplyr::left_join(validation_groups_df, by = "consensus_annotation") |>
    dplyr::left_join(gene_exp_df, by = c("barcodes", "library_id")) |>
    # add marker gene information (associated validation group annotation, gene observed count, percent tissues)
    # account for the same gene being present in multiple cell types
    dplyr::left_join(markers_df, by = "ensembl_gene_id", relationship = "many-to-many")

  # prep for plots
  # get total number of cells per final annotation group
  total_cells_df <- consensus_df |>
    dplyr::select(library_id, barcodes, broad_celltype_group) |>
    dplyr::distinct() |>
    dplyr::count(broad_celltype_group, name = "total_cells")

  # table with one row per unique broad cell type/ marker gene combination
  # first all cells in with the same broad_celltype_group (determined based on consensus_annotation) are grouped together
  # then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
  # do this for every possible marker gene/ validation group assignment
  # second we calculate the mean expression and mean percentage of all marker genes in a given validation group (this value is used only in the second section of the report)
  group_stats_df <- consensus_df |>
    # for each assigned cell type/marker gene combo get total detected and mean expression
    # group by both broad group and validation group to account for genes that are expressed in more than one cell type
    dplyr::summarize(
      .by = c("broad_celltype_group", "ensembl_gene_id", "validation_group_annotation"),
      detected_count = sum(detected),
      mean_exp = mean(logcounts)
    ) |>
    # add in validation group for marker genes
    # this includes all possible marker genes and all possible validation group assignments
    dplyr::left_join(markers_df, by = c("ensembl_gene_id", "validation_group_annotation"), relationship = "many-to-many") |> # add total cells
    dplyr::left_join(total_cells_df, by = c("broad_celltype_group")) |>
    # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild
    dplyr::filter(total_cells > 50) |>
    dplyr::mutate(
      # get total percent expressed
      percent_exp = (detected_count / total_cells) * 100,
      # account for NA/unknowns and set axes order
      broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown") |>
        factor(levels = c(unique(markers_df$validation_group_annotation), "unknown"))
    )

  # get list of celltypes to keep and assign colors
  celltype_groups <- group_stats_df |>
    dplyr::pull(broad_celltype_group) |>
    unique() |>
    as.character()

  # filter markers to those that are actually relevant
  # we will only plot the marker genes for cell types that are part of the assigned broad validation group for this group of samples
  # we don't care about plotting marker genes for cell types that aren't present here
  filtered_markers_df <- markers_df |>
    dplyr::filter(
      validation_group_annotation %in% celltype_groups,
      gene_symbol %in% group_stats_df$gene_symbol
    )

  # specify x axis order for dotplot
  marker_gene_order <- filtered_markers_df |>
    dplyr::pull(gene_symbol)

  # set order for cell types
  celltype_order <- unique(filtered_markers_df$validation_group_annotation)

  # filter out low expressed genes
  dotplot_df <- group_stats_df |>
    dplyr::filter(mean_exp > 0, percent_exp > 10) |>
    dplyr::arrange(broad_celltype_group) |>
    # add a label for the plot
    dplyr::mutate(y_label = as.factor(glue::glue("{broad_celltype_group} ({total_cells})"))) |>
    # remove marker genes that aren't present in final annotations and set x axis order
    dplyr::filter(gene_symbol %in% marker_gene_order) |>
    dplyr::mutate(
      # set orders of gene symbol and validation groups
      y_label = factor(y_label, levels = unique(y_label)),
      gene_symbol = factor(gene_symbol, levels = marker_gene_order),
      validation_group_annotation = factor(validation_group_annotation, levels = celltype_order)
    )


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

  suppressMessages(duckplyr::methods_restore()) # back to dplyr outside the function
  return(combined_plot)
}
