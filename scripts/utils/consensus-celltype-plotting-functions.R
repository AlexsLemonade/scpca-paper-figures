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
#' @param dotplot_size_range Vector specifying the point size range for dotplots, with a default of `c(1,6)` which matches the `ggplot2` default
#'
#' @returns Dot plot with summarized expression of marker genes for consensus cell types
marker_gene_dotplot <- function(
    sample_ids,
    consensus_results_dir,
    validation_groups_df,
    markers_df,
    celltype_colors, 
    dotplot_size_range = c(1,6) 
){
  
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

  # we do not want to show groups without markers in the plot
  # include NA_character_ here to ensure we keep the unknown cells, since we can't directly filter on that value in duckplyr
  allowed_groups <- c(unique(markers_df$validation_group_annotation), NA_character_)

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
  
  # check which broad cell type groups do not have marker genes
  # only consider cell types with > 50 cells
  all_diagnosis_celltype_groups <- total_cells_df |> 
    dplyr::filter(total_cells > 50) |>
    dplyr::pull(broad_celltype_group) |>
    unique()
  
  missing_groups <- setdiff(all_diagnosis_celltype_groups, allowed_groups) |> 
    paste(collapse = ",")
  
  # print out a message to show the cell types that are not visible
  # add these cell types to the legends
  if(missing_groups != ""){
    message(glue::glue("The following cell types will not be shown in this plot: {missing_groups}")) 
  } else {
    message("All cell types will be shown")
  }
  

  # table with one row per unique broad cell type/ marker gene combination
  # first all cells in with the same broad_celltype_group (determined based on consensus_annotation) are grouped together
  # then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
  # do this for every possible marker gene/ validation group assignment
  # second we calculate the mean expression and mean percentage of all marker genes in a given validation group (this value is used only in the second section of the report)
  group_stats_df <- consensus_df |>
    # only keep groups for which we have marker genes
    dplyr::filter(broad_celltype_group %in% allowed_groups) |>
    # for each assigned cell type/marker gene combo get total detected and mean expression
    # group by both broad group and validation group to account for genes that are expressed in more than one cell type
    dplyr::summarize(
      .by = c("broad_celltype_group", "ensembl_gene_id", "validation_group_annotation"),
      detected_count = sum(detected),
      mean_exp = mean(logcounts)
    ) |>
    # add in validation group for marker genes
    # this includes all possible marker genes and all possible validation group assignments
    dplyr::left_join(markers_df, by = c("ensembl_gene_id", "validation_group_annotation"), relationship = "many-to-many") |> 
    # add total cells
    dplyr::left_join(total_cells_df, by = c("broad_celltype_group")) |>
    # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild
    dplyr::filter(total_cells > 50) |>
    dplyr::mutate(
      # get total percent expressed
      percent_exp = (detected_count / total_cells) * 100,
      # account for NA/unknowns and set axes order
      broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown") |>
        factor(levels = unique(c(names(celltype_colors), "unknown"))) |>
        # ensure unknown is at the end
        forcats::fct_relevel("unknown", after = Inf)
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
    ) |> 
    # ensure order matches the order of the legend
    dplyr::mutate(
      validation_group_annotation = factor(validation_group_annotation, levels = names(celltype_colors))
    ) |> 
    dplyr::arrange(validation_group_annotation)

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
    dplyr::mutate(y_label = glue::glue("{broad_celltype_group} ({total_cells})")) |>
    # remove marker genes that aren't present in final annotations and set x axis order
    dplyr::filter(gene_symbol %in% marker_gene_order) |>
    dplyr::mutate(
      # set orders of gene symbol and validation groups
      y_label = factor(y_label, levels = rev(unique(y_label))),
      gene_symbol = factor(gene_symbol, levels = marker_gene_order),
      validation_group_annotation = factor(validation_group_annotation, levels = celltype_order)
    )

  combined_plot <- make_dotplot(
    dotplot_df,
    cell_type_label = "Broad cell type annotation",
    celltype_colors,
    dotplot_size_range
  )
  suppressMessages(duckplyr::methods_restore()) # back to dplyr outside the function
  return(combined_plot)
}


#' Helper function for creating the dot plots directly from the pre-formatted dataframes
#'
#' @param dotplot_df Data frame with y_label, gene_symbol, mean_exp, percent_exp, 
#'   and validation_group_annotation columns to use for plotting
#' @param cell_type_label Label to use for y axis of dotplot
#' @param celltype_colors Colors to use for cell type annotation bar
#' @param dotplot_size_range Vector specifying the point size range for dotplots, with a default of `c(1,6)` which matches the `ggplot2` default
#'
#' @returns dotplot with annotation bar
make_dotplot <- function(
    dotplot_df,
    cell_type_label,
    celltype_colors,
    dotplot_size_range = c(1,6)
) {
  
  # make dotplot with marker gene exp
  dotplot <- ggplot(dotplot_df, aes(y = y_label, x = gene_symbol, color = mean_exp, size = percent_exp)) +
    geom_point() +
    scale_size(range = dotplot_size_range) + 
    scale_color_viridis_c(option = "magma") +
    facet_grid(cols = vars(validation_group_annotation), scales = "free", space = "free") +
    guides(
      color = guide_colorbar(order = 1), 
      size = guide_legend(order = 2) 
    ) +
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
      y = cell_type_label,
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

#' Prep data frame to use for creating stacked bar plots showing all cell types 
#'
#' @param celltype_files List of files containing consensus cell type results
#' @param validation_groups_df Data frame with validation group assignments 
#'
#' @returns
#' @export
#'
#' @examples
create_celltype_summary <- function(
    celltype_files,
    validation_groups_df
){
  
  # read in consensus files and create data frame
  consensus_df <- celltype_files |> 
    purrr::map(readr::read_tsv) |> 
    dplyr::bind_rows()
  
  # get celltype summary for stacked bar chart 
  # need to add in validation groups here and do summary by validation group 
  consensus_df <- consensus_df |>
    # add in broad cell type group which is used for plotting
    # groups similar cell types together
    dplyr::left_join(validation_groups_df, by = "consensus_annotation") |> 
    # remove any PDX samples
    dplyr::filter(sample_type == "patient tissue") |> 
    # add in unknown for plotting 
    dplyr::mutate(broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown"))
  
  # get total cell count and number of assigned cell types per library
  totals_df <- consensus_df |> 
    dplyr::group_by(library_id) |> 
    dplyr::summarize(
      total_cells_per_library = dplyr::n()
    ) 
  
  # get summary stats for each cell type in each library  
  summary_df <- consensus_df |> 
    dplyr::left_join(totals_df, by = "library_id") |> 
    dplyr::group_by(project_id, library_id, sample_id, broad_celltype_group) |> 
    dplyr::summarize(total_cells_per_annotation = dplyr::n(),
                     total_cells_per_library = unique(total_cells_per_library),
                     percent_cells_annotation = round((total_cells_per_annotation / total_cells_per_library) * 100, 2)) |>
    dplyr::ungroup()
  
  # order by total % of annotated cells 
  # get a vector of library ids ordered by total percentage annotated
  library_levels <- summary_df |> 
    dplyr::filter(broad_celltype_group != "unknown") |> 
    dplyr::group_by(library_id) |> 
    dplyr::summarize(
      total_percent_annotated = sum(total_cells_per_annotation)/unique(total_cells_per_library)
    ) |>
    dplyr::arrange(desc(total_percent_annotated)) |> 
    dplyr::pull(library_id)
  
  # reorder by total percentage annotated 
  summary_df <- summary_df |> 
    dplyr::mutate(
      library_id = forcats::fct_relevel(library_id, library_levels),
      broad_celltype_group = forcats::fct_relevel(broad_celltype_group, "unknown", after = Inf) |> 
        forcats::fct_rev()
    ) |>
    unique()  
  
  return(summary_df)
}


#' Prep data frame to use for creating stacked bar plots showing immune cell types
#'   with an emphasis on T and myeloid cell types 
#'
#' @param celltype_files List of files containing consensus cell type results
#' @param all_immune_celltypes Vector of all consensus immune cell types possible
#' @param tcell_celltypes Vector of T cell types to plot 
#' @param myeloid_celltypes Vector of myeloid cell types to plot 
#' @param frac_immune_threshold Threshold of fraction of immune cells in library required to include a library in the figure
#'
#' @returns Summarized data frame for input to plotting
create_immune_celltype_summary <- function(
    celltype_files,
    all_immune_celltypes, 
    tcell_celltypes, 
    myeloid_celltypes, 
    frac_immune_threshold
){
  
  # read in consensus files and create data frame
  consensus_df <- celltype_files |> 
    purrr::map(readr::read_tsv) |> 
    dplyr::bind_rows()
  
  # Determine libraries to remove as those with < frac_immune_threshold
  #  fraction of immune cells out of total
  remove_libraries <- consensus_df |>
    # Keep only the immune cells and remove PDX
    dplyr::filter(
      sample_type == "patient tissue"
    ) |>
    dplyr::mutate(
      category = ifelse(
        consensus_annotation %in% c(tcell_celltypes, myeloid_celltypes), 
        "immune", 
        "other"
      )
    )  |>
    dplyr::count(library_id, category) |>
    tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0) |>
    dplyr::rowwise() |>
    dplyr::mutate(frac_immune = immune / (other + immune)) |>
    dplyr::filter(frac_immune < frac_immune_threshold) |>
    dplyr::pull(library_id)
  
  # subset to only immune celltypes, and add column for plotting
  consensus_df <- consensus_df |>
    # Keep only the immune cells and remove PDX
    dplyr::filter(
      consensus_annotation %in% all_immune_celltypes, 
      sample_type == "patient tissue"
    ) |>
    # Create immune_celltype_group label with value "other" when cells are not T or myeloid types
    dplyr::mutate(
      immune_celltype_group = ifelse(
        consensus_annotation %in% c(tcell_celltypes, myeloid_celltypes), 
        consensus_annotation, 
        "other"
      )
    ) 
  
  # get total cell count and number of assigned cell types per library
  totals_df <- consensus_df |> 
    dplyr::group_by(library_id) |> 
    dplyr::summarize(
      total_cells_per_library = dplyr::n()
    ) 
  
  # get summary stats for each cell type in each library  
  summary_df <- consensus_df |> 
    dplyr::left_join(totals_df, by = "library_id") |> 
    # remove libraries with insufficient cells
    dplyr::filter(!(library_id %in% remove_libraries)) |>
    dplyr::group_by(project_id, library_id, sample_id, immune_celltype_group) |> 
    dplyr::summarize(total_cells_per_annotation = dplyr::n(),
                     total_cells_per_library = unique(total_cells_per_library),
                     percent_cells_annotation = round((total_cells_per_annotation / total_cells_per_library) * 100, 2)) |>
    dplyr::ungroup()
  
  
  # Determine the order for immune cell categories based on:
  # - myeloid and t-cell types should be grouped together
  # - within each group, cell types should be ordered based on _overall frequency_
  # - finally, "other" should be first
  immune_factor_order <- summary_df |>
    dplyr::filter(immune_celltype_group != "other") |>
    # add up all the fractions as a proxy for overall frequency
    dplyr::group_by(immune_celltype_group) |>
    dplyr::summarize(total_frac = sum(percent_cells_annotation)) |>
    # assign groupings so we can order by them
    dplyr::mutate(
      immune_group = ifelse(immune_celltype_group %in% tcell_celltypes, "tcell", "myeloid")
    ) |>
    dplyr::group_by(immune_group) |>
    dplyr::arrange(desc(total_frac), .by_group = TRUE) |>
    dplyr::pull(immune_celltype_group)
  immune_factor_order <- c("other", immune_factor_order)
  
  # order by % of myeloid cells 
  # get a vector of library ids ordered by total percentage annotated
  library_levels <- summary_df |> 
    dplyr::filter(immune_celltype_group %in% myeloid_celltypes) |> 
    dplyr::group_by(library_id) |> 
    dplyr::summarize(
      myeloid_frac = sum(total_cells_per_annotation)/unique(total_cells_per_library)
    ) |>
    dplyr::arrange(desc(myeloid_frac)) |> 
    dplyr::pull(library_id)
  
  # reorder by total percentage annotated 
  summary_df <- summary_df |> 
    dplyr::mutate(
      library_id = forcats::fct_relevel(library_id, library_levels),
      immune_celltype_group = forcats::fct_relevel(immune_celltype_group, immune_factor_order)
    ) |>
    unique()  
  
  return(summary_df)
}


#' Stacked bar chart showing the percentage of cells annotated as each annotation, with optional faceting
#' Each column is a library ID and the fill of each bar corresponds to the percent of that sample annotated as that cell type
#'
#' @param df Data frame to use for plotting. Must have `fill_column`, `facet_variable` (if used), `library_id`, and `percent_cells_annotation` as columns
#' @param fill_column Column to use for determing fill color of each bar
#' @param celltype_colors Named vector of cell types and colors, names should match values in `fill_column`
#' @param fill_label Label for fill column to show on the legend
#' @param y_label Label for y-axis
#' @param lumped_label Label for cells that are shown in grey and should be last in the legend order (e.g., unknown or other)
#' @param facet_variable Column to use for faceting, default is NULL 
#' @param x_axis_text_size Size of x-axis text, default is 4
#' @param facet_col Number of columns to use in faceting, default is 2
#' @param facet_row Number of rows to use in faceting, default is 1
#' @param legend_position Where to put the legend, default is "right"
#'
#' @returns
#' @export
#'
#' @examples
stacked_barchart <- function(
    df, 
    fill_column,
    celltype_colors, # named vector where names match the values in fill_column
    fill_label = "Broad cell type annotation", 
    y_label = "Percent of cells",
    lumped_label = "unknown",
    facet_variable = NULL, # use for faceting HGG vs. LGG 
    x_axis_text_size = 4, 
    facet_col = 2,
    facet_row = 1,
    legend_position = "right"
){
  
  # make sure colors are named properly 
  stopifnot(
    "Names of celltype_colors must match values in fill_column" = all(df[[fill_column]] %in% names(celltype_colors))
  )
  
  plot_breaks <- c(
    setdiff(
      levels(df[[fill_column]]), 
      lumped_label
    ), 
    lumped_label
  )
  
  # if size is 0, we want element_blank()
  if (x_axis_text_size == 0) {
    axis_text_x_theme <- element_blank()
  } else {
    axis_text_x_theme <- element_text(angle = 60, hjust = 1, vjust = 1, size = x_axis_text_size)
  }
  
  barchart <- ggplot(df) + 
    aes(
      x = library_id, 
      y = percent_cells_annotation, 
      fill = !!sym(fill_column)
    ) +
    geom_col() + 
    scale_y_continuous(expand = c(0,0)) +
    # make sure unknown is last but all other legend order is based on when it appears 
    scale_fill_manual(
      values = celltype_colors,
      breaks = plot_breaks
    ) +
    theme_classic() +
    theme(
          axis.text.x = axis_text_x_theme,
          strip.background = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
          # add a square around each of the plots
          panel.background = element_rect(colour = "black", linewidth=0.5),
          axis.title = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 12),
          legend.position = legend_position) +
    labs(
      x = "", 
      y = y_label,
      fill= fill_label
    )
  
  if(!is.null(facet_variable)){
    barchart <- barchart +
      facet_wrap(vars(!!sym(facet_variable)), scales = "free_x", ncol = facet_col)
  }
  
  return(barchart)
  
}

#' Bar chart faceted by diagnosis group
#' Facets are dynamic based on the number of libraries in a diagnosis group
#'
#' @param plot_df Data frame to use for plotting. Must have `fill_column`, `facet_variable` (if used), `library_id`, and `percent_cells_annotation` as columns
#' @param fill_column Column to use for determing fill color of each bar
#' @param celltype_colors Named vector of cell types and colors, names should match values in `fill_column`
#' @param fill_label Label for fill column to show on the legend
#' @param lumped_label Label for cells that are shown in grey and should be last in the legend order (e.g., unknown or other)
#' @param diagnosis_column Column to use for faceting, default is "diagnosis_lumped"
#'
#' @returns
#' @export
#'
#' @examples
diagnosis_group_barchart <- function(
    plot_df, 
    fill_column, 
    celltype_colors, # named vector where names match the values in fill_column
    fill_label = "Broad cell type annotation", 
    lumped_label = "unknown",
    diagnosis_column = "diagnosis_lumped"
  ){
  
  plot_df <- plot_df |> 
    dplyr::group_by(!!sym(diagnosis_column)) |> 
    # add the number of libraries for each diagnosis
    dplyr::mutate(
      unique_library_count = dplyr::n_distinct(library_id)
    ) |> 
    dplyr::ungroup() |> 
    # set number of columns for faceted plot based on library number
    dplyr::mutate(
      ncol = dplyr::case_when(
        # numbers chosen to optimize aesthetics across all groups
        unique_library_count > 45 ~ 1,
        unique_library_count > 29 ~ 2,
        .default = 3
      )
    )
  
  # create a plot just to extract the legend
  legend_plot <- ggplot(plot_df) +
    aes(
      x = library_id, 
      y = percent_cells_annotation, 
      fill = !!sym(fill_column)
    ) +
    geom_col() +
    scale_fill_manual(
      values = celltype_colors,
      breaks = c(setdiff(unique(plot_df[[fill_column]]), lumped_label), lumped_label)
    ) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
      fill = fill_label
    )
  
  # extract the legend and make sure its a ggplot object 
  shared_legend <- ggpubr::get_legend(legend_plot) |> 
    ggplotify::as.ggplot()
  
  # split dataframe by number of columns needed in faceting
  split_df <- split(plot_df, plot_df$ncol, drop = TRUE)
  
  # create a list of plots, one for each number of columns in faceting 
  barchart_list <- split_df |> 
    purrr::map(\(df){
      
      faceted_plot <- stacked_barchart(
        df, 
        fill_column = fill_column, 
        celltype_colors = celltype_colors, 
        fill_label = fill_label,
        facet_variable = diagnosis_column, 
        facet_col = unique(df$ncol),
        x_axis_text_size = 0, # blank x-axis text
        legend_position = "none" # don't include legends 
      )
      
      return(faceted_plot)
      
    })
  
  # combine into one plot and add the legend 
  combined_barchart <- patchwork::wrap_plots(barchart_list, ncol = 1)  / shared_legend
  
  return(combined_barchart)
  
} 
