# This script generates a ridge plot of CNV events across cell types for a Neuroblastoma sample
# Code is adapted from the cell type QC report:
# https://github.com/AlexsLemonade/scpca-nf/blob/f8ee2a65887eaeb3ad137b64291984c3342e27f0/templates/qc_report/infercnv_qc.rmd#L129

# load project
renv::load()

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
})
theme_set(theme_bw())

# Set up -----------------------------------------------------------------------

project_id <- "SCPCP000004"
sample_id <- "SCPCS000112"
library_id <- "SCPCL000130"
data_dir <- here::here("s3_files", project_id, sample_id)

sce_file <- file.path(data_dir, glue::glue("{library_id}_processed.rds"))
cnv_palette_file <- here::here("palettes", "nb-cnv-palette.tsv")

output_pdf <- here::here("figures", "pdfs", "Fig5D_nb-ridgeplot.pdf")

# Prepare data for plotting ----------------------
sce <- readr::read_rds(sce_file)

celltype_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::select(infercnv_total_cnv, celltype = consensus_celltype_annotation) |>
  dplyr::mutate(
    celltype = stringr::str_to_lower(celltype), # be consistent with other paper figures
    celltype = forcats::fct_lump_n(
      celltype, 7, other_level = "all remaining cell types", ties.method = "first"
    ),
    celltype = forcats::fct_infreq(celltype),
    celltype = forcats::fct_relevel(
      celltype, "unknown", "all remaining cell types", after = Inf
    ), 
    celltype = forcats::fct_rev(celltype)
  )

# Make and export plot -----------------------------

nb_ridgeplot <- ggplot(celltype_df) +
  aes(x = infercnv_total_cnv, y = celltype, fill = celltype) +
  ggridges::geom_density_ridges2(scale = 0.9) +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  labs(
    x = "Total CNV events per cell",
    y = "Consensus cell type annotation"
  ) +
  theme(legend.position = "none")

ggsave(output_pdf, nb_ridgeplot, width = 4.5, height = 4.5)
