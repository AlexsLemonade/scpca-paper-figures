# This script is used to generate the heatmap comparing submitter annotations to SingleR/ CellAssign

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)

# Set heatmap padding option
ComplexHeatmap::ht_opt(TITLE_PADDING = grid::unit(0.2, "in"))

# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCS000251")

# define file paths to downloaded files 
processed_sce_file <- file.path(local_results_dir, "SCPCL000498_processed.rds")

# read in sce objects 
processed_sce <- readr::read_rds(processed_sce_file)

# define output file paths 
plots_dir <- here::here("figures", "pngs") 
output_plot_file <- file.path(plots_dir, "FigS4C_submitter-heatmap.png")

# source in helper functions for plotting
function_file <- here::here("scripts", "utils", "celltype-plot-helper-functions.R")
source(function_file)

# Prep data for plotting -------------------------------------------------------
# Create data frame of cell types
celltype_df <- create_celltype_df(processed_sce)

available_celltypes <- c("SingleR", "CellAssign")

jaccard_submitter_matrices <- available_celltypes |>
  stringr::str_to_lower() |>
  purrr::map(\(name) {
    make_jaccard_matrix(
      celltype_df,
      "submitter_celltype_annotation",
      glue::glue("{name}_celltype_annotation")
    )
  }) |>
  purrr::set_names(available_celltypes)

# Plot -------------------------------------------------------------------------

# set up png to save heatmap
png(output_plot_file, width = 9, height = 10, units = 'in', res = 300)

# create concatenated heatmaps 
jaccard_submitter_matrices |>
  purrr::imap(
    \(celltype_mat, celltype_method) {
      ComplexHeatmap::Heatmap(
        t(celltype_mat), # transpose because matrix rows are in common & we want a vertical arrangement
        col = circlize::colorRamp2(c(0, 1), colors = c("white", "darkslateblue")),
        border = TRUE, 
        ## Row parameters
        cluster_rows = FALSE,
        row_title = celltype_method,
        row_title_gp = grid::gpar(fontsize = 12),
        row_title_side = "right",
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = 10),
        ## Column parameters
        cluster_columns = FALSE,
        column_title = "Submitter-provided annotations",
        column_title_gp = grid::gpar(fontsize = 12),
        column_names_side = "bottom",
        column_names_gp = grid::gpar(fontsize = 10),
        column_names_rot = 90,
        ## Legend parameters
        heatmap_legend_param = list(
          title = "Jaccard index",
          direction = "vertical",
          legend_width = unit(1.5, "in")
        ),
        show_heatmap_legend = celltype_method == "SingleR",
      )
    }) |>
  # concatenate vertically into HeatmapList object
  purrr::reduce(ComplexHeatmap::`%v%`) |> 
  ComplexHeatmap::draw(
    heatmap_legend_side = "right",
    # add a margin to the heatmap so labels don't get cut off
    padding = unit(c(2, 20, 2, 2), "mm")
  )

dev.off()

