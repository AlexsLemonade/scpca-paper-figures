#!/usr/bin/env Rscript

# This file contains helper functions for plotting consensus cell types 

# packages required for functions
library(ggplot2)
library(patchwork)
library(data.table)


save_broad_group_dts <- function(
    celltype_file,
    gene_exp_file,
    validation_groups_dt,
    markers_dt,
    scratch_dir
){
  
  # make a folder for each library ID
  library_id <- stringr::str_remove(basename(celltype_file), "_processed_consensus-cell-types.tsv.gz$")
  output_dir <- file.path(scratch_dir, library_id)
  fs::dir_create(output_dir)
  
  # read in gene exp
  gene_dt <- fread(gene_exp_file)
  
  # read in cell types and join with validation groups and gene expression
  celltype_dt <- fread(celltype_file) |> 
    dplyr::left_join(validation_groups_dt, by = "consensus_annotation") |> 
    dplyr::left_join(gene_dt, by = c("barcodes", "library_id")) |>
    # add marker gene information (associated validation group annotation, gene observed count, percent tissues)
    # account for the same gene being present in multiple cell types
    dplyr::left_join(markers_dt, by = "ensembl_gene_id", relationship = "many-to-many") |>
    dplyr::mutate(detected = logcounts > 0,
                  # make sure we don't lose the NA when splitting
                  broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown")) 
  
  # split by celltype group 
  split_dt <- celltype_dt |> 
    split(celltype_dt$broad_celltype_group)
  
  # save individual files for each cell type group for a given library 
  names(split_dt) |> 
    purrr::walk(\(broad_group){
      group_name <- stringr::str_replace_all(broad_group, " ", "-")
      file_path <- file.path(output_dir, glue::glue("{group_name}_combined.tsv"))
      fwrite(split_dt[[broad_group]], file_path, sep = "\t", quote = FALSE, na = "NA")
      return(file_path)
    })
}

make_group_stats_dt <- function(
      broad_group,
      markers_dt, 
      scratch_dir
) {
  
  # first replace spaces and define output file 
  broad_group <- stringr::str_replace_all(broad_group, " ", "-")
  output_file <- file.path(scratch_dir, glue::glue("{broad_group}_stats.tsv"))
  
  broad_group_files <- list.files(scratch_dir, pattern = glue::glue("{broad_group}_combined.tsv"), recursive = TRUE, full.names = TRUE)
  
  if(length(broad_group_files) > 0 ){
    
    consensus_dt <- broad_group_files |> 
      purrr::map(fread) |> 
      data.table::rbindlist(fill = TRUE, use.names = TRUE) 
    
    message(glue::glue("Read in celltype files for {broad_group}"))
    
    # prep for plots 
    # get total number of cells per final annotation group 
    total_cells <- consensus_dt |> 
      dplyr::select(library_id, barcodes, broad_celltype_group) |> 
      dplyr::distinct() |> 
      dplyr::pull(barcodes) |> 
      length()
    
    # table with one row per unique broad cell type/ marker gene combination 
    # first all cells in with the same broad_celltype_group (determined based on consensus_annotation) are grouped together
    # then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
    # do this for every possible marker gene/ validation group assignment 
    # second we calculate the mean expression and mean percentage of all marker genes in a given validation group (this value is used only in the second section of the report)
    group_stats_df <- consensus_dt |> 
      # for each assigned cell type/marker gene combo get total detected and mean expression
      # group by both broad group and validation group to account for genes that are expressed in more than one cell type
      dplyr::group_by(broad_celltype_group, ensembl_gene_id, gene_symbol, validation_group_annotation) |>
      dplyr::summarize(
        detected_count = sum(detected),
        mean_exp = mean(logcounts),
        .groups = "drop"
      ) |> 
      dplyr::mutate(total_cells = total_cells) |> 
      # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild 
      dplyr::filter(total_cells > 50) |> 
      dplyr::rowwise() |> 
      dplyr::mutate(
        # get total percent expressed
        percent_exp = (detected_count/total_cells) * 100
      ) 
    
    readr::write_tsv(group_stats_df, output_file)
    
    message(glue::glue("Exported stats files for {broad_group}"))
    
    return(output_file)
    
  } 
  
}


#' Dot plot showing expression of marker genes across assigned cell types
#'
#' @param sample_ids Character vector of ScPCA sample ids to include in plot
#' @param consensus_results_dir Directory where results from cell-type-consensus module of OpenScPCA-nf lives
#' @param validation_groups_dt Data frame assigning consensus cell types to broader validation groups
#' @param markers_dt Data frame with marker genes for each cell type
#' @param celltype_colors Named vector of colors to use for each broader validation group
#'
#' @returns Dot plot with summarized expression of marker genes for consensus cell types
marker_gene_dotplot <- function(
    sample_ids,
    scratch_dir, 
    consensus_results_dir,
    validation_groups_dt,
    markers_dt,
    celltype_colors
){
  
  # # list all cell type assignments files
  # consensus_results_files <- list.files(consensus_results_dir, pattern = "_processed_consensus-cell-types\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE) 
  # celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]
  # 
  # # gene expression used for dot plots 
  # gene_exp_files <- list.files(consensus_results_dir, pattern = "_processed_marker-gene-expression\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  # gene_exp_files <- gene_exp_files[basename(dirname(gene_exp_files)) %in% sample_ids]
  # 
  # # read in files 
  # consensus_dt <- celltype_files |>
  #   purrr::map(fread) |>
  #   data.table::rbindlist(fill = TRUE, use.names = TRUE) 
  # 
  # gene_exp_dt <- gene_exp_files |>
  #   purrr::map(fread) |>
  #   data.table::rbindlist(fill = TRUE, use.names = TRUE)
  # 
  # message("Consensus cell type results read in")
  # 
  # # Join all consensus results and marker gene info
  # consensus_dt <- consensus_dt |>
  #   # add in broad cell type group which is used for plotting
  #   # groups similar cell types together
  #   dplyr::left_join(validation_groups_dt, by = "consensus_annotation") |>
  #   dplyr::left_join(gene_exp_dt, by = c("barcodes", "library_id")) |>
  #   # add marker gene information (associated validation group annotation, gene observed count, percent tissues)
  #   # account for the same gene being present in multiple cell types
  #   dplyr::left_join(markers_dt, by = "ensembl_gene_id", relationship = "many-to-many") |>
  #   dplyr::mutate(detected = logcounts > 0)
  # 
  # # remove extra data frame
  # rm(gene_exp_dt)
  # gc()
  # 
  # # prep for plots 
  # # get total number of cells per final annotation group 
  # total_cells_df <- consensus_dt |> 
  #   dplyr::select(library_id, barcodes, broad_celltype_group) |> 
  #   dplyr::distinct() |> 
  #   dplyr::count(broad_celltype_group, name = "total_cells") 
  # 
  # # table with one row per unique broad cell type/ marker gene combination 
  # # first all cells in with the same broad_celltype_group (determined based on consensus_annotation) are grouped together
  # # then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
  # # do this for every possible marker gene/ validation group assignment 
  # # second we calculate the mean expression and mean percentage of all marker genes in a given validation group (this value is used only in the second section of the report)
  # group_stats_df <- consensus_dt |> 
  #   # for each assigned cell type/marker gene combo get total detected and mean expression
  #   # group by both broad group and validation group to account for genes that are expressed in more than one cell type
  #   dplyr::group_by(broad_celltype_group, ensembl_gene_id, validation_group_annotation) |>
  #   dplyr::summarize(
  #     detected_count = sum(detected),
  #     mean_exp = mean(logcounts)
  #   ) |> 
  #   # add in validation group for marker genes
  #   # this includes all possible marker genes and all possible validation group assignments 
  #   dplyr::left_join(markers_dt, by = c("ensembl_gene_id", "validation_group_annotation"), relationship = "many-to-many") |>
  #   # now get the mean expression/ mean percentage across all marker genes for a given validation group
  #   # here the broad_celltype_group is the final assigned annotation for that group of cells 
  #   # the validation_group_annotation refers to the cell type that marker gene is associated with 
  #   dplyr::group_by(broad_celltype_group, validation_group_annotation) |> 
  #   dplyr::mutate(
  #     # calculate mean expression/detected across all markers for a specific group 
  #     all_markers_mean_exp = mean(mean_exp),
  #     all_markers_detected_count = mean(detected_count)
  #   ) |>  # add total cells
  #   dplyr::left_join(total_cells_df, by = c("broad_celltype_group")) |> 
  #   # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild 
  #   dplyr::filter(total_cells > 50) |> 
  #   dplyr::rowwise() |> 
  #   dplyr::mutate(
  #     # get total percent expressed
  #     percent_exp = (detected_count/total_cells) * 100,
  #     all_markers_percent_exp = (all_markers_detected_count/total_cells) * 100, 
  #     # account for NA/unknowns and set axes order
  #     broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown") |> 
  #       factor(levels =  c(unique(markers_dt$validation_group_annotation), "unknown"))
  #   ) 
  # 
  # # no longer need this and it takes up space 
  # rm(consensus_dt)
  # gc()
  
  # list all cell type assignments files
  consensus_results_files <- list.files(consensus_results_dir, pattern = "_processed_consensus-cell-types\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE) 
  celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]
  
  # gene expression used for dot plots 
  gene_exp_files <- list.files(consensus_results_dir, pattern = "_processed_marker-gene-expression\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  gene_exp_files <- gene_exp_files[basename(dirname(gene_exp_files)) %in% sample_ids]
  
  # combine cell type, gene exp, validation groups, and marker gene lists and save individual files to save on memory 
  purrr::walk2(celltype_files, gene_exp_files, 
               \(type_file, gene_file){
                 save_broad_group_dts(type_file, gene_file, validation_groups_dt, markers_dt, scratch_dir)
               })
  
  # get a list of all validation groups 
  broad_groups <- c(unique(markers_dt$validation_group_annotation), "unknown")
  
  all_stats_files <- broad_groups |> 
    purrr::map(\(group){
      make_group_stats_dt(group, markers_dt, scratch_dir)
    }) |> 
    unlist()
  
  group_stats_df <- all_stats_files |> 
    purrr::map(fread) |> 
    data.table::rbindlist(fill = TRUE, use.names = TRUE) |> 
    dplyr::mutate(
      # account for NA/unknowns and set axes order
      broad_celltype_group = factor(broad_celltype_group, levels = broad_groups)
    )
  
  # get list of celltypes to keep and assign colors 
  celltype_groups <- group_stats_df |> 
    dplyr::pull(broad_celltype_group) |> 
    unique() |>
    as.character()
  
  # filter markers to those that are actually relevant 
  # we will only plot the marker genes for cell types that are part of the assigned broad validation group for this group of samples
  # we don't care about plotting marker genes for cell types that aren't present here 
  filtered_markers_df <- markers_dt |> 
    dplyr::filter(validation_group_annotation %in% celltype_groups,
                  gene_symbol %in% group_stats_df$gene_symbol)

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

  return(combined_plot)
}
