# Generates figures and manuscript numbers tables for bulk analysis

renv::load()

# load any libraries 
library(ggplot2)

# Set up -----------------------------------------------------------------------

analysis_dir <- here::here("analysis", "pseudobulk-bulk-prediction")
result_dir <- file.path(analysis_dir, "results")

tables_dir <- here::here("manuscript-numbers")
bulk_numbers_tsv <- file.path(tables_dir, "bulk-analysis-counts.tsv")

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

sample_metadata_file <- here::here("s3_files", "scpca-library-metadata.tsv")


# First, create the manuscript-numbers table ----------------------------------------------------------------

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


# Export -----------------------------------------------------------------------
readr::write_tsv(bulk_table, bulk_numbers_tsv)
