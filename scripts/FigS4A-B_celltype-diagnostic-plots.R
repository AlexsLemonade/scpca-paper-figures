# This script is used to generate example diagnostic plots for CellAssign and SingleR 
# Much of this code comes from the supplemental cell type report: 
# https://github.com/AlexsLemonade/scpca-nf/blob/main/templates/qc_report/celltypes_supplemental_report.rmd

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12)
    )
)

# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCS000001")

# define file paths to downloaded files 
processed_sce_file <- file.path(local_results_dir, "SCPCL000001_processed.rds")

# read in sce objects 
processed_sce <- readr::read_rds(processed_sce_file)

# define output file paths 
pdf_dir <- here::here("figures", "pdfs") 
singler_diagnostic_pdf_file <- file.path(pdf_dir, "FigS4A_singler-diagnostic.pdf")
cellassign_diagnostic_pdf_file <- file.path(pdf_dir, "FigS4B_cellassign-diagnostic.pdf")

# source in helper functions for plotting
function_file <- here::here("scripts", "utils", "celltype-plot-helper-functions.R")
source(function_file)

# Create data frame of cell types
celltype_df <- create_celltype_df(processed_sce)

# Prep SingleR data for plotting -------------------------------------------------------

# extract scores into matrix
singler_scores <- metadata(processed_sce)$singler_result$scores

# Create data frame for plotting with delta median and the full *non-pruned* cell labels
delta_median_df <- tibble::tibble(
  delta_median = rowMaxs(singler_scores) - rowMedians(singler_scores),
  # Need to grab the non-pruned label for this plot
  full_labels = metadata(processed_sce)$singler_result$labels,
  # if pruned.labels are NA ==> low confidence
  # so, negate for this variable:
  confident = !is.na(metadata(processed_sce)$singler_result$pruned.labels)
) |>
  dplyr::mutate(
    confident = ifelse(confident, "High-quality", "Low-quality")
  )


# we use inner_join b/c celltype_df does NOT contain "Unknown cell type", which
#  we do not want to display here
delta_median_df <- delta_median_df |>
  dplyr::inner_join(
    tibble::tibble(
      full_labels = celltype_df$singler_celltype_ontology,
      celltype = celltype_df$singler_celltype_annotation
    ) |> dplyr::distinct()
  ) |>
  dplyr::select(-full_labels)

# add column with ordered levels with wrapped labels for visualization
delta_median_df$annotation_wrapped <- factor(
  delta_median_df$celltype,
  # rev() so large groups are at the TOP of the plot
  levels = rev(levels(delta_median_df$celltype)),
  labels = rev(stringr::str_wrap(levels(delta_median_df$celltype), 30))
)

# Subset the data to just confident points for median+/-IQR
delta_median_confident_df <- delta_median_df |> 
  dplyr::filter(confident == "High-quality")

# SingleR diagnostic plot ------------------------------------------------------

# Plot delta_median across celltypes colored by pruning
singler_diagnostic_plot <- ggplot(delta_median_df) +
  aes(
    x = delta_median,
    y = annotation_wrapped,
    shape = confident,
    alpha = confident
  ) +
  ggforce::geom_sina(
    size = 0.8,
    color = "black", # will get applied to all confident points and non-confident outline
    fill = "white", # will apply to non-confident fill only
    position = position_dodge(width = 0.05) # Keep both types of points mostly in line
  ) +
  # Handle points aesthetics:
  #  confident are closed black with alpha = 0.5
  #  not confident are open black with alpha = 1
  scale_shape_manual(values = c(19, 21)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(
    x = "Delta median statistic",
    y = "Cell type annotation",
    shape = "Cell type annotation quality"
  ) +
  # add median diamond for confident points only
  stat_summary(
    data = delta_median_confident_df,
    color = "red",
    geom = "point",
    fun = "median",
    shape = 18,
    size = 2.25,
    alpha = 0.9
  ) +
  guides(
    alpha = "none",
    shape = guide_legend(override.aes = list(size = 1.5, alpha = 0.55))
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(singler_diagnostic_pdf_file, singler_diagnostic_plot, height = 7, width = 7)

# CellAssign plotting ----------------------------------------------------------

# define bandwidth for all calculations
density_bw <- 0.03

# find the maximum density across all distributions, and
#  save the maximum for determining geom_segment height
y_max <- celltype_df$cellassign_max_prediction |>
  split(celltype_df$cellassign_celltype_annotation) |>
  # make sure we get rid of any small groups
  purrr::discard(\(x) sum(is.finite(x)) <= 2) |>
  # remove any NA's that may have slipped in
  purrr::map_dbl(
    \(x) max(density(x, bw = density_bw, na.rm = TRUE)$y)
  ) |>
  max(na.rm = TRUE)

# add count to celltype_df for setting alpha and yend values
celltype_df <- celltype_df |>
  dplyr::add_count(cellassign_celltype_annotation)

# make the plot!
cellassign_diagnostic_plot <- ggplot(celltype_df) +
  aes(x = cellassign_max_prediction) +
  geom_density(
    bw = density_bw,
    fill = "grey65",
    linewidth = 0.25,
    bounds = c(0, 1)
  ) +
  geom_segment(
    aes(
      # set alpha to vary based on the number of points in the row such that
      #  rows with more points are more transparent
      alpha = pmax(0.2, 1 - 0.01 * n),
      xend = cellassign_max_prediction,
      # set yend as either 0 for rows with many points, or y_max/2.5 for
      #  rows with few points
      yend = ifelse(n > 5, 0, y_max / 2.5),
      y = -Inf
    ),
    color = "blue"
  ) +
  labs(
    x = "Probability of annotated cell type",
    y = "Cell type annotation"
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_alpha_identity() +
  facet_grid(
    rows = vars(cellassign_celltype_annotation),
    switch = "y"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 1
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.02, "in")
  )

ggsave(cellassign_diagnostic_pdf_file, cellassign_diagnostic_plot, width = 7, height = 7)
