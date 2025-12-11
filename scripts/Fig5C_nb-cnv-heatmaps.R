# This script generates heatmaps of CNV events for a Neuroblastoma sample

# load project
renv::load()

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ComplexHeatmap)
  library(ggplot2)
})

# Set up -----------------------------------------------------------------------

project_id <- "SCPCP000004"
sample_id <- "SCPCS000112"
library_id <- "SCPCL000130"
data_dir <- here::here("s3_files", project_id, sample_id)

sce_file <- file.path(data_dir, glue::glue("{library_id}_processed.rds"))
cnv_palette_file <- here::here("palettes", "nb-cnv-palette.tsv")

# define output files; we need a separate file for each heatmap since 
# each is already vertically concatenated, and ComplexHeatmap can't horizontally concatenate them
output_chr1_file <- here::here("figures", "pdfs", "Fig5C_chr1-heatmap.pdf")
output_chr11_file <- here::here("figures", "pdfs", "Fig5C_chr11-heatmap.pdf")
output_chr17_file <- here::here("figures", "pdfs", "Fig5C_chr17-heatmap.pdf")
legend_file <- here::here("figures", "pdfs", "Fig5C_legend.pdf")

# chromosomes to make heatmaps for
chrs_to_plot <- c(1, 11, 17) |>
  purrr::set_names() |>
  purrr::map(
    \(chr) paste0(rep(chr, each = 2), c("p","q"))
  )

# define heatmap pdfs as a list for use within purrr::map
output_files <- list(
  output_chr1_file, 
  output_chr11_file, 
  output_chr17_file
) |>
  purrr::set_names(names(chrs_to_plot))

# Define helper functions ----------------

#' Create a data frame of CNV events per cell, for a given CNV type
#'
#' @param metadata_table infercnv metadata table with HMM results
#' @param celltypes_df data frame with OpenScPCA cell types
#' @param arms_to_keep vector of chromosome arms to keep
#' @param cnv_type which type of CNV to tabulate, either "dupli" or "loss"
#'
#' @returns Data frame with CNV events per cell, where loss CNVs are negative
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
    # need rowwise for _only_ this operation
    dplyr::rowwise() |>
    # record losses as negative for plotting
    dplyr::mutate(cnv = ifelse(cnv_type == "loss", -1*cnv, cnv))
}


#' Make a heatmap of CNV events for a given category of cells
#'
#' @param df Data frame with per cell cnv events and cell type categories
#' @param category Cell type category to plot
#'
#' @returns a ComplexHeatmap object
make_heatmap <- function(df, category) {
  
  # create matrix for plotting heatmap
  heatmap_mat <- df |>
    dplyr::filter(cell_category == category) |>            
    dplyr::select(-cell_category) |>
    tidyr::pivot_wider(
      names_from = chr,
      values_from = cnv
    ) |>
    dplyr::arrange(barcodes) |>
    tibble::column_to_rownames("barcodes") |>
    as.matrix()
  
  # height should be scaled based on the number of cells
  heatmap_height <- nrow(heatmap_mat) * 0.01
  
  # make the heatmap
  Heatmap(
    heatmap_mat,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE, 
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 10),
    column_names_centered = TRUE,
    show_heatmap_legend = FALSE, # we're making our own legend
    border = TRUE, # helps to see chr with few CNV events
    height = grid::unit(heatmap_height, "mm")
  )
}

# Prepare data for plotting ----------------------
sce <- readr::read_rds(sce_file)

# define palette
cnv_palette <- readr::read_tsv(cnv_palette_file) |>
  dplyr::mutate(
    value = dplyr::case_when(
      cnv_type == "gain" ~ 1,
      cnv_type == "loss" ~ -1,
      cnv_type == "none" ~ 0
    )
  ) |>
  dplyr::arrange(value)
col_fun <- circlize::colorRamp2(
  cnv_palette$value,
  cnv_palette$color
)

# define cell type categories
unknown_celltypes <- c("Unknown", "openscpca-excluded")
celltypes_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::select(barcodes, openscpca_celltype_annotation) |>
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
event_df <- dplyr::bind_rows(dupli_df, loss_df) 

# Set order of cells:
# each cell will be represented with two rows: gain event, and then loss event
# this code orders barcodes in order gain/loss
barcode_levels <- paste0(rep(colnames(sce), each = 2), c("_dupli", "_loss"))
event_df <- event_df |>
  dplyr::mutate(
    barcodes = factor(barcodes, levels = barcode_levels)
  ) 

# Make and export a heatmap per chromosome ------------------------------
heatmap_list <- chrs_to_plot |>
  purrr::imap(
  \(chr_arms, chr) {

    chr_cnv_df <- event_df |>
      dplyr::filter(chr %in% chr_arms)
    
    mal_heatmap <- make_heatmap(chr_cnv_df, "malignant")
    normal_heatmap <- make_heatmap(chr_cnv_df, "normal")
    unknown_heatmap <-  make_heatmap(chr_cnv_df, "unknown")
    
    stacked <- mal_heatmap %v% normal_heatmap %v% unknown_heatmap
    
    pdf(output_files[[as.character(chr)]], height = 4.5, width = 0.75)
    draw(stacked, ht_gap = grid::unit(c(2, 4), "mm")) # bigger gap before unknown
    dev.off()
  }
) 

# Make a legend ------------------------------------------

# dummy data frame for getting a line segment
plot_3col  <- data.frame(
  cnv_type = rep(cnv_palette$cnv_type, each = 2),
  x = 1:6,
  y = 6:11
) |>
  ggplot() +
  aes(x, y, color = cnv_type) +
  geom_line(linewidth = 1) +
  theme_void() +
  scale_color_manual(
    name = "CNV type",
    values = cnv_palette$color |> rlang::set_names(cnv_palette$cnv_type)) +
  theme(legend.key = element_rect(fill = "grey75", color = NA))
legend_3col <- cowplot::get_legend(plot_3col)
ggsave(legend_file, legend_3col, width = 0.75, height = 1)
