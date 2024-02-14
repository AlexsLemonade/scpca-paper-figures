# This script is used to generate a heatmap comparing SingleR to CellAssign cell type annotations

# load project
renv::load()

library(SingleCellExperiment)

# Set heatmap padding option
ComplexHeatmap::ht_opt(TITLE_PADDING = grid::unit(0.2, "in"))


# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCS000001")

# define file paths to downloaded files 
processed_sce_file <- file.path(local_results_dir, "SCPCL000001_processed.rds")

# read in sce objects 
processed_sce <- readr::read_rds(processed_sce_file)

# define output file paths 
figure_dir <- here::here("figures", "pdfs")
output_file <- file.path(figure_dir, "Fig4B_singler-cellassign-heatmap.pdf")

# source in helper functions for plotting
function_file <- here::here("scripts", "utils", "celltype-plot-helper-functions.R")
source(function_file)

# Prep data for plotting -------------------------------------------------------
# Create data frame of cell types
celltype_df <- create_celltype_df(processed_sce)
  
# Calculate jaccard matrix
singler_cellassign_matrix <- make_jaccard_matrix(
  celltype_df,
  "singler_celltype_annotation",
  "cellassign_celltype_annotation"
)

# Plot -------------------------------------------------------------------------

# heatmap comparing singleR to cellassign annotations
heatmap <- ComplexHeatmap::Heatmap(
  # transpose because matrix rows are in common & we want a vertical arrangement
  t(singler_cellassign_matrix),
  col = circlize::colorRamp2(c(0, 1), colors = c("white", "darkslateblue")),
  border = TRUE,
  ## Row parameters
  cluster_rows = FALSE,
  row_title = "CellAssign annotations",
  row_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "right",
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 10),
  ## Column parameters
  cluster_columns = FALSE,
  column_title = "SingleR annotations",
  column_title_gp = grid::gpar(fontsize = 12),
  column_names_side = "bottom",
  column_names_gp = grid::gpar(fontsize = 10),
  # ensure column labels fit in PDF export
  column_names_max_height = grid::unit(8, "cm"),
  ## Legend parameters
  heatmap_legend_param = list(
    title = "Jaccard index",
    direction = "vertical",
    legend_width = grid::unit(1.5, "in")
  )
) |>
ComplexHeatmap::draw(
  heatmap_legend_side = "right"
)

# save heatmap to pdf
pdf(output_file, width = 9, height = 9, useDingbats = FALSE)
ComplexHeatmap::draw(heatmap)
dev.off()
