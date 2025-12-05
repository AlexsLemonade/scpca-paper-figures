# Let's make a shit ton of heat maps
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ComplexHeatmap)
})


outdir <- "so-many-heatmaps/contenders/query-reference"
fs::dir_create(outdir)

chr_levels <- paste0(rep(1:22, each = 2), c("p","q"))
col_fun <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("dodgerblue4", "white", "firebrick")
)

chrs_to_test <- 1:22
keep_chrs <- paste0(rep(chrs_to_test, each = 2), c("p","q"))

get_percell_cnv_df <- function(cnv_df, col_starts){
  cnv_df |>
    tidyr::pivot_longer(
      starts_with(col_starts),
      names_to = "chr",
      values_to = "cell_cnv"
    ) |>
    dplyr::mutate(
      chr = stringr::str_replace(chr, col_starts, ""), 
      chr = factor(chr, levels = chr_levels)
    ) |>
    dplyr::select(library_id, cell_id, chr, cell_cnv)
}

prep_for_heatmap <- function(df, category, chrs, cell_id_levels) {
  df |>
    dplyr::filter(
      cell_category == category,
      chr %in% chrs
    ) |>            
    dplyr::select(-cell_category) |>
    tidyr::pivot_wider(
      names_from = chr,
      values_from = cell_cnv
    ) |>
    dplyr::mutate(cell_id = factor(cell_id, levels = cell_id_levels)) |>
    dplyr::arrange(cell_id)|>
    tibble::column_to_rownames("cell_id") |>
    as.matrix()
}

make_heatmap <- function(mat, col_fun = col_fun) {
   Heatmap(
    mat,
    name = "CNV",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE, 
   # column_title = title,
    height = grid::unit(nrow(mat) *0.005, "mm")
  )
}


sce <- readRDS("~/Desktop/SCPCP000004_merged.rds")
sce_metadata <- metadata(sce)

sce_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::select(cell_id, library_id, is_infercnv_reference, openscpca_celltype_annotation) |>
  dplyr::mutate(
    #cell_category = dplyr::case_when(
    #  openscpca_celltype_annotation == "Neuroendocrine" ~ "malignant", 
    #  openscpca_celltype_annotation %in% c("Unknown", "openscpca-excluded") ~ "Unknown", 
    #  .default = "non-malignant")
    cell_category = ifelse(
      is_infercnv_reference, 
      "non-malignant", 
      "malignant")
    ) |>
  dplyr::filter(cell_category!="Unknown") |>
  tibble::as_tibble()


processed_libraries <- sce_df |>
  dplyr::filter(is_infercnv_reference) |>
  dplyr::count(library_id) |>
  dplyr::filter(n > 100) |>
  dplyr::pull(library_id)


cnv_df <- processed_libraries |>
  purrr::set_names() |>
  purrr::map(
    \(library_id) {
      sce_metadata$library_metadata[[library_id]]$infercnv_table
    }
  ) |>
  purrr::list_rbind(names_to = "library_id") |>
  dplyr::mutate(cell_id = glue::glue("{library_id}-{barcodes}")) |>
  dplyr::select(-barcodes)




gain_df <- get_percell_cnv_df(cnv_df, "has_dupli_chr") |>
  dplyr::left_join(sce_df, by = c("library_id", "cell_id")) |>
  dplyr::mutate(cell_id = glue::glue("{cell_id}_gain"))
loss_df <- get_percell_cnv_df(cnv_df, "has_loss_chr") |>
  dplyr::left_join(sce_df, by = c("library_id", "cell_id")) |>
  dplyr::mutate(cell_id = glue::glue("{cell_id}_loss"),
                cell_cnv = -cell_cnv)
event_df <- dplyr::bind_rows(gain_df, loss_df) |>
  dplyr::filter(chr %in% keep_chrs)

c("SCPCL000130") |> 
  purrr::map(
    \(lib_id) {
      library_cnv_df <- event_df |>
        dplyr::filter(library_id == lib_id) |>
        dplyr::select(cell_id, cell_category, chr, cell_cnv)
      
      all_cell_ids <- library_cnv_df$cell_id |>
        stringr::str_replace("_\\w{4}$", "") |>
        unique()
      
      cell_id_levels <- paste0(rep(all_cell_ids, each = 2), c("_gain", "_loss"))
      
      out_height <- nrow(library_cnv_df) / 100000

      #### map over chromosomes ######
      heatmap_list <- purrr::map(
        chrs_to_test, \(chr) {
          keep_me <- c(
            glue::glue("{chr}p"), 
            glue::glue("{chr}q")
          )
          
          mal_heatmap <- library_cnv_df |>
            prep_for_heatmap("malignant", keep_me, cell_id_levels) |>
            make_heatmap()
          

          normal_heatmap <- library_cnv_df |>
            prep_for_heatmap("non-malignant", keep_me, cell_id_levels) |>
            make_heatmap()
          
          
          stacked <- mal_heatmap %v% normal_heatmap
          fs::dir_create(glue::glue("{outdir}/{lib_id}"))
          pdf(glue::glue("{outdir}/{lib_id}/heatmap_{chr}.pdf"), width = 1, height = 3)
          draw(stacked, ht_gap = grid::unit(5, "mm"))
          dev.off()
          
        }
      ) 
    }
  )
# 
# 
# Neuroblastoma	1	p	blue	!!!!!!!!!!	
# Neuroblastoma	1	q	red	
# Neuroblastoma	2	p	red	
# Neuroblastoma	3	p	blue	
# Neuroblastoma	4	p	blue	
# Neuroblastoma	6	q	blue	
# Neuroblastoma	7	q	red	
# Neuroblastoma	8	p	blue	
# Neuroblastoma	11	q	blue	!!!!!!!!!!!!!!!!
# Neuroblastoma	12	q	red	
# Neuroblastoma	14	q	blue	
# Neuroblastoma	17	q	red	!!!!!!!!!!!!!!!
# Neuroblastoma	19	p	blue	
# Neuroblastoma	19	q	blue	
# Neuroblastoma	22	q	blue	