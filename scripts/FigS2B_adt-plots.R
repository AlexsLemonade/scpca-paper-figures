# This script is used to create the mini versions of the plots included in the ADT report
# we use the results from SCPCS000216/SCPCL000290

renv::load()

library(ggplot2)
library(SingleCellExperiment)

theme_set(
  theme_classic() +
    theme(
      strip.background = element_rect(fill = "transparent"),
      # no axis ticks or labels
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # add a square around each of the plots
      panel.background = element_rect(colour = "black", linewidth=0.5),
      aspect.ratio = 1
    )
)

set.seed(2024)


# Set up -----------------------------------------------------------------------

# define file paths and read in sce objects
data_dir <- here::here("s3_files", "SCPCS000216")

filtered_sce_file <- file.path(data_dir, "SCPCL000290_filtered.rds")
filtered_sce <- readr::read_rds(filtered_sce_file)

processed_sce_file <- file.path(data_dir, "SCPCL000290_processed.rds")
processed_sce <- readr::read_rds(processed_sce_file)

output_plot_file <- here::here("figures", "pngs", "FigS2B_adt-plots.png")

# Filtered ADT plot ------------------------------------------------------------

# grab coldata from filterd object 
filtered_coldata_df <- colData(filtered_sce) |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode")

filter_levels <- c("Keep", "Filter (RNA & ADT)", "Filter (RNA only)", "Filter (ADT only)")
filtered_coldata_df <- filtered_coldata_df |>
  # add column to represent filtering on both RNA and ADT counts
  dplyr::mutate(filter_summary = dplyr::case_when(
    scpca_filter == "Keep" & adt_scpca_filter == "Keep" ~ filter_levels[1],
    scpca_filter == "Remove" & adt_scpca_filter == "Remove" ~ filter_levels[2],
    scpca_filter == "Remove" & adt_scpca_filter == "Keep" ~ filter_levels[3],
    scpca_filter == "Keep" & adt_scpca_filter == "Remove" ~ filter_levels[4],
    TRUE ~ "Error"
  ),
  # make sure Keep is first 
  filter_summary = forcats::fct_relevel(filter_summary, "Keep"),
  keep_column = ifelse(
    stringr::str_starts(filter_summary, "Keep"),
    "Keep",
    "Remove"
  ))

# faceted plot to show cells belonging to each filter group
filtered_plot <- ggplot(filtered_coldata_df, aes(x = detected, y = subsets_mito_percent, color = keep_column)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = "Number of genes detected",
    y = "Mitochondrial percentage"
  ) +
  facet_wrap(vars(filter_summary)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Top ADTs ---------------------------------------------------------------------

# Calculate variance for each ADT
adt_var <- altExp(processed_sce) |>
  logcounts() |>
  apply(1, var, na.rm = TRUE)

# Get the top 4
top_adts <- adt_var[order(adt_var, decreasing = TRUE)[1:4]] |>
  names()

# ADT expression density -------------------------------------------------------

# grab expression for top ADTs from counts
var_adt_exp_df <- logcounts(altExp(processed_sce))[top_adts, ] |>
  # as.matrix needs to come _first_
  as.matrix() |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode") |>
  # combine all ADTs into a single column for  faceting
  tidyr::pivot_longer(
    !barcode,
    names_to = "ADT",
    values_to = "adt_expression"
  )

adt_density_plot <- ggplot(var_adt_exp_df, aes(x = adt_expression, fill = ADT)) +
  geom_density() +
  facet_wrap(vars(ADT), nrow = 2) +
  labs(x = "Log-normalized ADT expression",
       y = "Density") +
  theme(legend.position = "none")

# ADT expression UMAP ----------------------------------------------------------

# create data frame of UMAPs and expression
umap_df <- scuttle::makePerCellDF(processed_sce) |>
  tibble::rownames_to_column("barcode") |>
  dplyr::select(
    barcode,
    UMAP1 = UMAP.1,
    UMAP2 = UMAP.2
  ) |>
  # combine with gene expression
  dplyr::left_join(var_adt_exp_df, by = "barcode")

adt_umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = adt_expression)) +
  geom_point(alpha = 0.1, size = 0.1) +
  facet_wrap(vars(ADT)) +
  scale_color_viridis_c() +
  labs(
    color = "Log-normalized ADT expression"
  ) +
  # remove axis numbers and background grid
  scale_x_continuous(labels = NULL, breaks = NULL) +
  scale_y_continuous(labels = NULL, breaks = NULL) +
  coord_fixed() +
  theme(
    legend.position = "none"
  )

# Combine and export -----------------------------------------------------------

combined_plot <- patchwork::wrap_plots(list(filtered_plot, adt_density_plot, adt_umap_plot))

ggsave(output_plot_file, combined_plot)
