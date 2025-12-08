# This script generates heatmaps for Figure 5C

# load project
renv::load()

library(SingleCellExperiment)
library(ComplexHeatmap)


# Set up -----------------------------------------------------------------------

project_id <- "SCPCP000004"
sample_id <- "SCPCS000112"
library_id <- "SCPCL000130"
data_dir <- here::here("s3_files", project_id, sample_id)

sce_file <- file.path(data_dir, glue::glue("{library_id}_processed.rds"))
sce <- readr::read_rds(sce_file)

# define output files; we need a separate file for each heatmap since 
# each is already vertically concatenated, and ComplexHeatmap can't horizontally concatenate them
output_chr1_file <- here::here("figures", "pdfs", "Fig5C-chr1-heatmap.pdf")
output_chr11_file <- here::here("figures", "pdfs", "Fig5C-chr11-heatmap.pdf")
output_chr17_file <- here::here("figures", "pdfs", "Fig5C-chr17-heatmap.pdf")
legend_file <- here::here("figures", "pdfs", "Fig5C-legend.pdf")


gain_color <- "firebrick4"
loss_color <- "blue"
none_color <- "white"

col_fun <- circlize::colorRamp2(
  c(-1, 0, 1),
  c(loss_color, none_color, gain_color)
)


chrs_to_plot <- c(1, 11, 17) |>
  purrr::set_names() |>
  purrr::map(
    \(chr) paste0(rep(chr, each = 2), c("p","q"))
  )

# define as a list for use within purrr::map
output_files <- list(
  output_chr1_file, 
  output_chr11_file, 
  output_chr17_file
) |>
  purrr::set_names(names(chrs_to_plot))

# Define helper functions ----------------

# Used to create a data frame of either gain ("dupli") or loss events per cell, with cell types
prepare_percell_cnv_df <- function(metadata_table, celltypes_df, cnv_type) {
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
      barcode_type = glue::glue("{barcodes}_{cnv_type}")
    ) |>
    # need rowwise for _only_ this operation
    dplyr::rowwise() |>
    dplyr::mutate(cnv = ifelse(cnv_type == "loss", -1*cnv, cnv) # losses are recorded as negative
    )
}


make_heatmap <- function(df, category) {
  
  # create matrix for plotting heatmap
  heatmap_mat <- df |>
    dplyr::filter(cell_category == category) |>            
    dplyr::select(-cell_category) |>
    tidyr::pivot_wider(
      names_from = chr,
      values_from = cnv
    ) |>
    dplyr::arrange(barcode_type) |>
    tibble::column_to_rownames("barcode_type") |>
    as.matrix()
  
  # height should be scaled based on the number of cells
  heatmap_height <- nrow(heatmap_mat) * 0.01
  
  # make the heatmap
  Heatmap(
    heatmap_mat,
    col = col_fun,
    name = "cnv",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE, # will be removed in illustrator 
    show_heatmap_legend = FALSE, # we're making our own legend
    height = grid::unit(heatmap_height, "mm")
  )
}




# Prepare data for plotting ----------------------
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
dupli_df <- prepare_percell_cnv_df(
  metadata(sce)$infercnv_table, celltypes_df, "dupli"
)
loss_df <- prepare_percell_cnv_df(
  metadata(sce)$infercnv_table, celltypes_df, "loss"
)
event_df <- dplyr::bind_rows(dupli_df, loss_df) |>
  dplyr::filter(chr %in% purrr::reduce(chrs_to_plot, c))

# Set order of cells: gain above loss
all_barcodes <- event_df$barcodes |>
  stringr::str_replace("_\\w{4}$", "") |>
  unique()
barcode_levels <- paste0(rep(all_barcodes, each = 2), c("_dupli", "_loss"))
event_df <- event_df |>
  dplyr::mutate(
    barcode_type = factor(barcode_type, levels = barcode_levels)
  ) |>
  dplyr::select(-barcodes)

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



