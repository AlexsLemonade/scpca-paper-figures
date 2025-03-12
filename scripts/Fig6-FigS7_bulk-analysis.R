# Generates figures and manuscript numbers tables for bulk analysis

renv::load()

# load any libraries 
library(ggplot2)
theme_set(theme_bw())



# Set up -----------------------------------------------------------------------

# significance threshold for odds ratios
alpha <- 0.05 

# Define projects intended for main vs supplementary figures
main_projects <- c("SCPCP000001", "SCPCP000002", "SCPCP000009")
si_projects   <- c("SCPCP000006", "SCPCP000017")

# File paths
analysis_dir <- here::here("analysis", "pseudobulk-bulk-prediction")
result_dir <- file.path(analysis_dir, "results")
tables_dir <- here::here("manuscript-numbers")
sample_metadata_file <- here::here("s3_files", "scpca-library-metadata.tsv")

# This file contains all samples which were _initially_ considered before 
# low-quality samples were removed
map_file <- file.path(analysis_dir, "data", "bulk-library-sample-ids.tsv")

# The TSV files that contain the geneset lists, needed for manuscript numbers
geneset_files <- list.files(
  path = file.path(analysis_dir, "data"), 
  pattern = "_panglao-genesets\\.tsv$", 
  full.names = TRUE
) |>
  purrr::set_names(
    \(x) {stringr::str_split_i(basename(x), pattern = "_", i = 1)}
  )

# The TSV files containing odds ratio results from overrepresentation analysis
odds_files <- list.files(
  path = result_dir,
  pattern = "_ORA-odds-ratios\\.tsv$",
  full.names = TRUE
)|>
  purrr::set_names(
    \(x) {stringr::str_split_i(basename(x), pattern = "_", i = 1)}
  )


# The TSV files containing modeled values
data_files <- list.files(
  path = file.path(result_dir, "models"),
  pattern = "-data_threshold-0\\.tsv$",
  full.names = TRUE
)|>
  purrr::set_names(
    \(x) {stringr::str_split_i(basename(x), pattern = "_", i = 1)}
  )


# Output files
bulk_numbers_tsv <- file.path(tables_dir, "bulk-analysis-counts.tsv")


# Functions --------------------------------------------------------------------

# Function to plot the odds ratio barplot
plot_odds_ratios <- function(df) {
  ggplot(df) +
    aes(
      x = tidytext::reorder_within(geneset_cell_type, odds_ratio, within = project_facet),
      y = odds_ratio, 
      fill = -log10(p_adj_bh)
    ) +
    geom_col() +
    scale_fill_viridis_c() +
    tidytext::scale_x_reordered() + # gets the labels back to only geneset_cell_type
    facet_wrap(vars(project_facet), scales = "free_x", nrow = 1) +
    labs(
      x = "Gene set (cell type)", 
      y = "Odds ratio", 
      fill = "-Log10 adjusted p-value (BH)"
    ) +
    theme(
      legend.position = "bottom", 
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8)
    )
}



# First, create the manuscript-numbers table -----------------------------------

# Define vector of low-quality samples which were removed
# Source: https://github.com/AlexsLemonade/scpca-paper-figures/blob/aa53bb48f2a7b1c9312033df25c912f25eb08147/analysis/pseudobulk-bulk-prediction/notebooks/build-assess-models.Rmd#L41-L53
excluded_samples <- c(
  "SCPCS000003", 
  "SCPCS000030", 
  "SCPCS000040", 
  "SCPCS000178", 
  "SCPCS000191", 
  "SCPCS000197", 
  "SCPCS000203", 
  "SCPCS000608"
)

# Read in data
sample_metadata <- readr::read_tsv(sample_metadata_file)
bulk_ids_df <- readr::read_tsv(map_file) |>
  # add in the project id too
  dplyr::inner_join(
    dplyr::select(sample_metadata, scpca_project_id, scpca_sample_id, scpca_library_id)
  ) |>
  # add indicator if low-quality
  dplyr::mutate(low_quality = scpca_sample_id %in% excluded_samples)

# Count how many samples per project were originally identified for use
total_table <- bulk_ids_df |>
  dplyr::count(scpca_project_id, name = "n_samples_total")

# Count how many samples per project were removed due to quality
excluded_table <- bulk_ids_df |>
  dplyr::filter(low_quality) |>
  dplyr::count(scpca_project_id, name = "n_samples_low_quality")

# Count the number of genesets used for ORA for each project
geneset_count_df <- geneset_files |>
  purrr::map(readr::read_tsv) |>
  purrr::list_rbind(names_to = "scpca_project_id") |>
  dplyr::select(scpca_project_id, cell_type_name) |>
  dplyr::distinct() |>
  dplyr::count(scpca_project_id, name = "n_ora_genesets")

# Combine all the count tables
bulk_table <- total_table |>
  dplyr::left_join(excluded_table) |>
  # replace NA with 0
  tidyr::replace_na(list(n_samples_low_quality = 0)) |>
  # create n_samples_used column
  dplyr::rowwise() |>
  dplyr::mutate(n_samples_used = n_samples_total - n_samples_low_quality) |>
  dplyr::ungroup() |>
  # final join
  dplyr::left_join(geneset_count_df)

# Finally, add a row with totals for manuscript writing convenience
bulk_table <- bulk_table |>
  tibble::add_row(
    scpca_project_id      = "total", 
    n_samples_total       = sum(bulk_table$n_samples_total), 
    n_samples_low_quality = sum(bulk_table$n_samples_low_quality), 
    n_samples_used        = sum(bulk_table$n_samples_used), 
    n_ora_genesets        = NA_integer_ # the sum of this quantity isn't relevant
  )


# Odds ratio panels -----------------------------------------------------------------

# Read in odds ratio results and filter to significant only
odds_df <- odds_files |>
  purrr::map(readr::read_tsv) |>
  purrr::list_rbind(names_to = "scpca_project_id") |>
  # filter to significant at 0.05
  dplyr::filter(p_adj_bh <= alpha)

# Prepare data frame for plotting  
odds_df <- odds_df |>
  # join with the number of samples in each project
  dplyr::inner_join(
    dplyr::select(bulk_table, scpca_project_id, n_samples_used)
  ) |>
  # create columns for plotting 
  dplyr::mutate(
    # faceting with the sample count info
    project_facet = glue::glue("{scpca_project_id} (N={n_samples_used})"), 
    # wrapped gene set names for space
    geneset_cell_type = stringr::str_wrap(geneset_cell_type, 25)
  )

# Main text panel:
odds_panel_main <- odds_df |>
  dplyr::filter(scpca_project_id %in% main_projects) |>
  plot_odds_ratios()


# SI panel:
odds_panel_si <- odds_df |>
  dplyr::filter(scpca_project_id %in% si_projects) |>
  plot_odds_ratios()




# Export -----------------------------------------------------------------------
readr::write_tsv(bulk_table, bulk_numbers_tsv)
