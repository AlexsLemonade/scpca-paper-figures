# This script generates UMAP panels 5A and 5B 

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)

# set default ggplot theme for UMAPs
theme_set(
  theme_classic() +
    theme(
      strip.background = element_rect(fill = "transparent", linewidth = 0.5),
      # no axis ticks or labels
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1
    )
)

# Set up -----------------------------------------------------------------------

project_id <- "SCPCP000004"
data_dir <- here::here("s3_files", project_id)

merged_sce_file <- file.path(data_dir, glue::glue("{project_id}_merged.rds"))
merged_sce_df <- readr::read_rds(merged_sce_file) |>
  # convert to data frame immediately for space
  scuttle::makePerCellDF(
    use.coldata = c("cell_id", "library_id", "infercnv_total_cnv", "openscpca_celltype_annotation"),
    use.dimred = "UMAP"
  ) |>
  # remove rownames
  tibble::as_tibble()

# define output files
output_panel_a_file <- here::here("figures", "pdfs", "Fig5A-umap-celltypes.pdf")
output_panel_b_file <- here::here("figures", "pdfs", "Fig5B-umap-infercnv.pdf")

  
# Prepare data for plotting ----------------------------------------------------
# Collapse cell types
merged_sce_df <- merged_sce_df |>
  dplyr::mutate(
    cell_category = dplyr::case_when(
      openscpca_celltype_annotation == "Neuroendocrine" ~ "Malignant cell", 
      openscpca_celltype_annotation %in% c("Unknown", "openscpca-excluded") ~ "Unknown",
      openscpca_celltype_annotation %in% c("Fibroblast", "Endothelial", "Schwann") ~ glue::glue("{openscpca_celltype_annotation} cell"),
      # everything else is a type of immune cell
      .default = "Immune cell"
    ),
    # determine factor order based on frequency, with Unknown last
    cell_category = forcats::fct_infreq(cell_category), 
    cell_category = forcats::fct_relevel(cell_category, "Unknown", after = Inf)
  )

# TODO this needs to be in its own file read in as input
cell_palette <- c(
  "Malignant cell" = "rosybrown4",
  "Immune cell" = "darkgreen",
  "Endothelial cell" = "#FC1CBF",
  "Fibroblast cell" = "#B10DA1", 
  "Schwann cell" = "chocolate1",
  "Unknown" = "gray70" # darker than validation palette since we have some standalone unknown clumps
)

# Panel A ------------------
celltype_umap <- ggplot(merged_sce_df) + 
  aes(x = UMAP.1, y = UMAP.2, color = cell_category) + 
  geom_point(alpha = 0.25, size = 0.25) + 
  scale_color_manual(values = cell_palette) + 
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 1.5))
  ) +
  labs(
    x = "UMAP1", 
    y = "UMAP2", 
    color = "OpenScPCA broad cell type"
  ) +
  theme(
    legend.position = "bottom", 
    legend.direction = "horizontal"
  ) 
  

# Panel B -----------------------------
infercnv_umap <- ggplot(merged_sce_df) + 
aes(x = UMAP.1, y = UMAP.2, color = infercnv_total_cnv) + 
  geom_point(alpha = 0.25, size = 0.25) + 
  scale_color_viridis_c(na.value = "grey70") +
  labs(
    x = "UMAP1", 
    y = "UMAP2", 
    color = "Total CNV"
  ) +
  theme(legend.position = "bottom") 


# Export --------------------------------
ggsave(output_panel_a_file, celltype_umap, width = 6, height = 6)
ggsave(output_panel_b_file, infercnv_umap, width = 6, height = 6)
