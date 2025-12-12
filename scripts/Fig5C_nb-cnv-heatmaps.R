# This script generates heatmaps of CNV events for a Neuroblastoma sample

# load project
renv::load()

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ComplexHeatmap)
})

# Set up -----------------------------------------------------------------------

project_id <- "SCPCP000004"
sample_id <- "SCPCS000112"
library_id <- "SCPCL000130"
data_dir <- here::here("s3_files", project_id, sample_id)

sce_file <- file.path(data_dir, glue::glue("{library_id}_processed.rds"))
cnv_palette_file <- here::here("palettes", "nb-cnv-palette.tsv")

output_heatmap_file <- here::here("figures", "pdfs", "Fig5C_cnv-heatmap.pdf")

# chromosomes to include in heatmap
chrs_to_plot <- paste0(rep(c(1, 11, 17), each = 2), c("p","q"))

# Define helper function ---------------------------------

#' Create a data frame of CNV events per cell, for a given CNV type
#'
#' @param metadata_table infercnv metadata table with HMM results
#' @param celltypes_df data frame with OpenScPCA cell types
#' @param arms_to_keep vector of chromosome arms to keep
#' @param cnv_type which type of CNV to tabulate, either "dupli" or "loss"
#'
#' @returns Data frame with CNV events per cell, where events are coded as discrete values
prepare_percell_cnv_df <- function(
    metadata_table, 
    celltypes_df, 
    arms_to_keep, 
    cnv_type
  ) {
  keep_col_starts <- glue::glue("has_{cnv_type}_chr")
  
  metadata_table |>
    dplyr::select(barcodes, starts_with(keep_col_starts)) |>
    tidyr::pivot_longer(
      -barcodes,
      names_to = "chr",
      values_to = "cnv"
    ) |>
    # combine with cell type information
    dplyr::left_join(celltypes_df, by = "barcodes") |>
    # prepare for plotting
    dplyr::mutate(
      chr = stringr::str_replace(chr, keep_col_starts, ""),
      barcodes = glue::glue("{barcodes}_{cnv_type}")
    ) |>
    # only keep relevant barcodes
    dplyr::filter(chr %in% arms_to_keep) |>
    # make it discrete
    dplyr::mutate(cnv = ifelse(cnv == 0, "none", cnv_type))
}


# Prepare data for plotting ----------------------
sce <- readr::read_rds(sce_file)

# define palette
cnv_palette <- readr::read_tsv(cnv_palette_file) |>
  tibble::deframe()

# define cell type categories
unknown_celltypes <- c("Unknown", "openscpca-excluded")
celltypes_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::mutate(
    cell_category = dplyr::case_when(
      openscpca_celltype_annotation == "Neuroendocrine" ~ "malignant", 
      openscpca_celltype_annotation %in% unknown_celltypes ~ "unknown",
      .default = "normal"
    )
  ) |>
  dplyr::select(barcodes, cell_category)


# Create data frame with per-cell events, with gain and loss on separate rows
all_arms <- purrr::reduce(chrs_to_plot, c)
dupli_df <- prepare_percell_cnv_df(
  metadata(sce)$infercnv_table, celltypes_df, all_arms, "dupli"
)
loss_df <- prepare_percell_cnv_df(
  metadata(sce)$infercnv_table, celltypes_df, all_arms, "loss"
)

# Set order of cells:
# each cell will be represented with two rows: gain event, and then loss event
# this code orders barcodes in order gain/loss
barcode_levels <- paste0(rep(colnames(sce), each = 2), c("_dupli", "_loss"))
event_df <- dupli_df |> 
  dplyr::bind_rows(loss_df) |>
  # pivot for making a matrix
  tidyr::pivot_wider(
    names_from = chr,
    values_from = cnv
  ) |>
  dplyr::mutate(
    barcodes = factor(barcodes, levels = barcode_levels),
    cell_category = factor(cell_category, levels = c("malignant", "normal", "unknown"))
  ) |>
  dplyr::arrange(cell_category, barcodes)

# Extract discrete matrix
event_mat <- event_df |>
  dplyr::select(-barcodes, -cell_category) |>
  as.matrix()

# Make and export the heatmap ----------------------------------

# scale height based on nubmer of cells
heatmap_height <- nrow(event_mat) * 0.01

# make the heatmap
ht <- ComplexHeatmap::Heatmap(
  event_mat,
  col = cnv_palette,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE, 
  column_names_rot = 0,
  column_names_gp = gpar(fontsize = 10),
  column_names_centered = TRUE,
  show_heatmap_legend = FALSE, # make separately
  border = TRUE, # helps to see chr with few CNV events
  # split by chr and cell type, but don't show the titles
  row_split    = event_df$cell_category,
  column_split = rep(c(1, 11, 17), each = 2),
  row_gap = unit(c(2, 4), "mm"), # more distance for unknown
  column_gap = unit(4, "mm"),
  column_title = NULL,
  row_title_gp = gpar(fontsize = 10),
  row_title_rot = 0,
  height = grid::unit(heatmap_height, "mm")
)

# make colors a bit transparent in the legend to look more similar to the plot itself
cnv_palette_alpha <- glue::glue("{cnv_palette}75")

# make the legend
lgd <- ComplexHeatmap::Legend(
  title = "CNV type",
  labels = c("Loss", "Neutral", "Gain"),
  legend_gp = gpar(fill = cnv_palette_alpha), 
  labels_gp = gpar(fontsize = 10),
  border = "gray10", 
  nrow = 1, 
  title_position = "leftcenter"
)

# Export
pdf(output_heatmap_file, height = 4.5, width = 4)
ComplexHeatmap::draw(
  ht,
  annotation_legend_list = list(lgd),
  merge_legend = TRUE, # needed because of splits
  heatmap_legend_side = "bottom"
)
dev.off()
