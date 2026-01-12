# This script generates UMAP panels of the Neuroblastoma project samples, colored by cell types and CNV events
# It also exports a table for manuscript-numbers giving the total number of labeled cells

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

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
palette_file <- here::here("palettes", "nb-annotation-palette.tsv")

merged_sce_file <- file.path(data_dir, glue::glue("{project_id}_merged.rds"))
merged_sce_df <- readr::read_rds(merged_sce_file) |>
  # convert to data frame immediately for space
  scuttle::makePerCellDF(
    use.coldata = c(
      "cell_id",
      "library_id",
      "infercnv_total_cnv",
      "openscpca_celltype_annotation",
      "consensus_celltype_annotation"
    ),
    use.dimred = "UMAP"
  ) |>
  # remove rownames
  tibble::as_tibble()

# define output file
# we'll export both panels at one for improved alignment
output_panel_ab_file <- here::here("figures", "pngs", "Fig5A-B_umap-celltypes-infercnv.png")
numbers_file <- here::here("manuscript-numbers", "SCPCP000004_labeled-cells.tsv")

# Prepare data for plotting ----------------------------------------------------
# Collapse cell types
merged_sce_df <- merged_sce_df |>
  dplyr::mutate(
    cell_category = dplyr::case_when(
      openscpca_celltype_annotation == "Neuroendocrine" ~ "malignant",
      openscpca_celltype_annotation %in%
        c("Unknown", "openscpca-excluded") ~ "unknown",
      openscpca_celltype_annotation %in%
        c(
          "Fibroblast",
          "Endothelial",
          "Schwann"
        ) ~ openscpca_celltype_annotation,
      # everything else is a type of immune cell
      .default = "immune cell"
    ),
    # tack on "cell" in a few spots
    cell_category = ifelse(
      cell_category %in% c("Endothelial", "Schwann", "malignant"),
      glue::glue("{cell_category} cell"),
      cell_category
    ),
    # make everything lowercase, consistent with other figures
    cell_category = stringr::str_to_lower(cell_category),
    # determine factor order based on frequency, with Unknown last
    cell_category = forcats::fct_infreq(cell_category),
    cell_category = forcats::fct_relevel(cell_category, "unknown", after = Inf)
  )

# prepare palette
cell_palette <- readr::read_tsv(palette_file) |> 
  tibble::deframe()


# Panel A ------------------
celltype_umap <- ggplot(merged_sce_df) +
  aes(x = UMAP.1, y = UMAP.2, color = cell_category) +
  geom_point(alpha = 0.25, size = 0.25) +
  scale_color_manual(values = cell_palette) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 3))
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "OpenScPCA broad\ncell type annotation"
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
combined_plot <- celltype_umap + infercnv_umap
ggsave(output_panel_ab_file, combined_plot, width = 10, height = 6)

# Calculate manuscript numbers and export
unknown_labels <- c("Unknown", "openscpca-excluded") # only Unknown for consensus, but same result
percent_openscpca <- sum(!(merged_sce_df$openscpca_celltype_annotation %in% unknown_labels)) / nrow(merged_sce_df)
percent_consensus <- sum(!(merged_sce_df$consensus_celltype_annotation %in% unknown_labels)) / nrow(merged_sce_df)

tibble::tribble(
  ~percent_labeled_openscpca, ~percent_labeled_consensus,
  percent_openscpca, percent_consensus
) |>
  readr::write_tsv(numbers_file)
