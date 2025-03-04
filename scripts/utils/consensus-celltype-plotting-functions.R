#!/usr/bin/env Rscript

# function to generat dot plots for a given set of samples 
marker_gene_dotplot <- function(
    sample_ids,
    consensus_results_dir,
    validation_groups_dt,
    markers_dt,
    celltype_colors
){
  
  # list consensus results 
  # cell type assignments 
  consensus_results_files <- list.files(consensus_results_dir, pattern = "_processed_consensus-cell-types\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE) 
  celltype_files <- consensus_results_files[basename(dirname(consensus_results_files)) %in% sample_ids]
  
  # gene expression used for dot plots 
  gene_exp_files <- list.files(consensus_results_dir, pattern = "_processed_marker-gene-expression\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  gene_exp_files <- gene_exp_files[basename(dirname(gene_exp_files)) %in% sample_ids]
  
  # read in files 
  consensus_dt <- celltype_files |> 
    purrr::map(fread) |>
    data.table::rbindlist(fill = TRUE, use.names = TRUE) 
  
  gene_exp_dt <- gene_exp_files |> 
    purrr::map(fread) |> 
    data.table::rbindlist(fill = TRUE, use.names = TRUE)
  
  consensus_dt <- consensus_dt |> 
    # add in broad cell type group which is used for plotting
    # groups similar cell types together 
    dplyr::left_join(validation_groups_dt, by = "consensus_annotation") |> 
    dplyr::left_join(gene_exp_dt, by = c("barcodes", "library_id")) |> 
    # add marker gene information (associated validation group annotation, gene observed count, percent tissues)
    # account for the same gene being present in multiple cell types 
    dplyr::left_join(markers_dt, by = "ensembl_gene_id", relationship = "many-to-many") |> 
    dplyr::mutate(detected = logcounts > 0)
  
  rm(gene_exp_dt)
  
  # prep for plots 
  # get total number of cells per final annotation group 
  total_cells_df <- consensus_dt |> 
    dplyr::select(library_id, barcodes, broad_celltype_group) |> 
    dplyr::distinct() |> 
    dplyr::count(broad_celltype_group, name = "total_cells") 
  
  # table with one row per unique broad cell type/ marker gene combination 
  # first all cells in with the same broad_celltype_group (determined based on consensus_annotation) are grouped together
  # then get the mean gene expression and total percentage of cells that express each marker gene across all cells in that group
  # do this for every possible marker gene/ validation group assignment 
  # second we calculate the mean expression and mean percentage of all marker genes in a given validation group (this value is used only in the second section of the report)
  group_stats_df <- consensus_dt |> 
    # for each assigned cell type/marker gene combo get total detected and mean expression
    # group by both broad group and validation group to account for genes that are expressed in more than one cell type
    dplyr::group_by(broad_celltype_group, ensembl_gene_id, validation_group_annotation) |>
    dplyr::summarize(
      detected_count = sum(detected),
      mean_exp = mean(logcounts)
    ) |> 
    # add in validation group for marker genes
    # this includes all possible marker genes and all possible validation group assignments 
    dplyr::left_join(markers_dt, by = c("ensembl_gene_id", "validation_group_annotation"), relationship = "many-to-many") |>
    # now get the mean expression/ mean percentage across all marker genes for a given validation group
    # here the broad_celltype_group is the final assigned annotation for that group of cells 
    # the validation_group_annotation refers to the cell type that marker gene is associated with 
    dplyr::group_by(broad_celltype_group, validation_group_annotation) |> 
    dplyr::mutate(
      # calculate mean expression/detected across all markers for a specific group 
      all_markers_mean_exp = mean(mean_exp),
      all_markers_detected_count = mean(detected_count)
    ) |>  # add total cells
    dplyr::left_join(total_cells_df, by = c("broad_celltype_group")) |> 
    # for plotting we're only going to look at any cell types with > 50 cells otherwise these plots can get wild 
    dplyr::filter(total_cells > 50) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      # get total percent expressed
      percent_exp = (detected_count/total_cells) * 100,
      all_markers_percent_exp = (all_markers_detected_count/total_cells) * 100, 
      # account for NA/unknowns and set axes order
      broad_celltype_group = tidyr::replace_na(broad_celltype_group, "unknown") |> 
        factor(levels =  c(unique(markers_dt$validation_group_annotation), "unknown"))
    ) 
  
  # no longer need this and it takes up space 
  rm(consensus_dt)
  
  # make a data frame that just has the unique genes 
  # unique genes are those that are only observed in one cell type based on observance in Cell Marker 
  # this was determined when creating the marker gene table from Cell Marker 
  unique_gene_df <- group_stats_df |> 
    # keep all 6 HPC genes 
    dplyr::filter(gene_observed_count == 1 | validation_group_annotation == "hematopoietic precursor cell")
  
  # get list of celltypes to keep and assign colors 
  celltype_groups <- group_stats_df |> 
    dplyr::pull(broad_celltype_group) |> 
    unique() |>
    as.character()
  
  # filter markers to those that are actually relevant 
  # we will only plot the marker genes for cell types that are part of the assigned broad validation group for this project
  # we don't care about plotting marker genes for cell types that aren't present here 
  # note that we will use this for both the dotplot and the heatmap
  filtered_markers_df <- markers_dt |> 
    dplyr::filter(validation_group_annotation %in% celltype_groups,
                  gene_symbol %in% group_stats_df$gene_symbol)
  
  # now only keep unique markers for the dotplot
  # except for hematopoietic precursor cells, where there are no unique markers
  unique_markers_df <- filtered_markers_df |> 
    dplyr::filter(gene_observed_count == 1 | validation_group_annotation == "hematopoietic precursor cell")
  
  # specify x axis order for dotplot
  marker_gene_order <- unique_markers_df |> 
    dplyr::pull(gene_symbol)
  
  # set order for cell types 
  celltype_order <- unique(unique_markers_df$validation_group_annotation)
  
  # filter out low expressed genes
  dotplot_df <- unique_gene_df |> 
    dplyr::filter(mean_exp > 0, percent_exp > 10) |> 
    dplyr::arrange(broad_celltype_group) |> 
    # add a label for the plot 
    dplyr::mutate(y_label = as.factor(glue::glue("{broad_celltype_group} ({total_cells})"))) |> 
    # remove marker genes that aren't present in final annotations and set x axis order 
    dplyr::filter(gene_symbol %in% marker_gene_order) |> 
    dplyr::mutate(
      # set orders of gene symbol and validation groups 
      gene_symbol = factor(gene_symbol, levels = marker_gene_order),
      validation_group_annotation = factor(validation_group_annotation, levels = celltype_order)
    )
  
  # make dotplot with marker gene exp
  dotplot <- ggplot(dotplot_df, aes(y = forcats::fct_rev(y_label), x = gene_symbol, color = mean_exp, size = percent_exp)) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      text = element_text(size = 14)
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
    scale_fill_manual(values = celltype_colors, breaks = levels(dotplot_df$validation_group_annotation)) +
    ggmap::theme_nothing() +
    theme(legend.position = "bottom") +
    labs(fill = "")
  
  combined_plot <- dotplot + color_bar +
    patchwork::plot_layout(ncol = 1, heights = c(4, 0.1)) 
  
  return(combined_plot)
}