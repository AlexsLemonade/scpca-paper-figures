# This script is used to generate panel A of the supplemental figure for CellAssign
# the output is a faceted UMAP showing the assigned cell types
# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)

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
      legend.key = element_blank()
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
umap_png_file <- file.path(png_dir, "FigS5A_cellassign-umap.png")

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
