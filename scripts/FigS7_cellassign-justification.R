# This script is used to generate the supplemental figure for CellAssign
# panel A is a faceted UMAP showing the assigned cell types
# panel B is a heatmap comparing CellAssign annotations to submitter annotations

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)

# Set heatmap padding option
ComplexHeatmap::ht_opt(TITLE_PADDING = grid::unit(0.2, "in"))

# Set default ggplot theme
theme_set(
  theme_classic() +
    theme(
      #plot.margin = margin(rep(20, 4)),
      strip.background = element_rect(fill = "transparent", linewidth = 0.5),
      # no axis ticks or labels
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # add a square around each of the plots
      panel.background = element_rect(colour = "black", linewidth = 0.5),
      aspect.ratio = 1,
      # remove boxes around legends 
      legend.key=element_blank()
    )
)

# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCS000264")

# define file paths to downloaded files
processed_sce_file <- file.path(local_results_dir, "SCPCL000490_processed.rds")

# read in sce objects
processed_sce <- readr::read_rds(processed_sce_file)

# define output file paths
# use png for UMAP 
png_dir <- here::here("figures", "pngs")
umap_png_file <- file.path(png_dir, "FigS7A_cellassign-umap.png")

# use pdf for heatmap 
pdf_dir <- here::here("figures", "pdfs")
heatmap_pdf_file <- file.path(pdf_dir, "FigS7B_cellassign-submitter-heatmap.pdf")

# source in helper functions for plotting
function_file <- here::here("scripts", "utils", "celltype-plot-helper-functions.R")
source(function_file)

# UMAP -------------------------------------------------------------------------
# Create data frame of cell types
celltype_df <- create_celltype_df(processed_sce)

# lump celltypes for cellassign 
celltype_df <- celltype_df |> 
  dplyr::mutate(cellassign_celltype_annotation_lumped = cellassign_celltype_annotation |> 
                  forcats::fct_lump_n(5, other_level = "All remaining cell types", ties.method = "first") |>
                  forcats::fct_infreq() |>
                  forcats::fct_relevel("Unknown cell type", "All remaining cell types", after = Inf)
  )

# faceted UMAP
faceted_umap <- ggplot(
  celltype_df,
  aes(x = UMAP1, y = UMAP2, color = cellassign_celltype_annotation_lumped)
) +
  # set points for all "other" points
  geom_point(
    data = dplyr::select(
      celltype_df, -cellassign_celltype_annotation_lumped
    ),
    color = "gray80",
    alpha = 0.5,
    size = 0.3
  ) +
  # set points for desired cell type
  geom_point(size = 0.3, alpha = 0.5) +
  facet_wrap(
    vars(cellassign_celltype_annotation_lumped),
    ncol = 3
  ) +
  scale_color_brewer(palette = "Dark2") +
  # remove axis numbers and background grid
  scale_x_continuous(labels = NULL, breaks = NULL) +
  scale_y_continuous(labels = NULL, breaks = NULL) +
  guides(
    color = guide_legend(
      title = "Cell types",
      # more visible points in legend
      override.aes = list(
        alpha = 1,
        size = 2
      )
    )
  ) +
  theme(legend.position = "bottom")

ggsave(umap_png_file, faceted_umap, width = 8, height = 8)

# Heatmap ----------------------------------------------------------------------

# get jaccard similarity index
jaccard_submitter_matrix <- make_jaccard_matrix(
  celltype_df,
  "submitter_celltype_annotation",
  "cellassign_celltype_annotation"
)

# heatmap comparing cellassign to submitter annotations
heatmap <- ComplexHeatmap::Heatmap(
  # transpose so submitter will be columns
  t(jaccard_submitter_matrix),
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
  column_title = "Submitter-provided annotations",
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
    heatmap_legend_side = "right", 
    # add a margin to the heatmap so labels don't get cut off
    padding = unit(c(2, 20, 2, 2), "mm")
  )

# save heatmap to pdf
pdf(heatmap_pdf_file, width = 9, height = 9, useDingbats = FALSE)
ComplexHeatmap::draw(heatmap)
dev.off()
