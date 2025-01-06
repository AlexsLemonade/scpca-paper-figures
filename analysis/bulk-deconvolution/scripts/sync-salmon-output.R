# This script is used to download data needed for bulk solid tumor deconvolution analysis
# All files will be saved in `../data/salmon-quant-files/<project id>/<sample id>/<library id>`
# 
# This script assumes that the ../../s3_files/scpca-<sample/library>-metadata.tsv files present. These files can be obtained with ../../scripts/figure_setup/sync-metadata.R.

renv::load()

# Define files and directories ---------
library_metadata_file <- here::here("s3_files", "scpca-library-metadata.tsv")
sample_metadata_file <- here::here("s3_files", "scpca-sample-metadata.tsv")
bulk_dir <- here::here("analysis", "bulk-deconvolution")
s3_path <- "s3://nextflow-ccdl-results/scpca-prod/checkpoints/salmon" #/library_id/quant.sf
data_path <- file.path(bulk_dir, "data", "salmon-quant-files")
fs::dir_create(data_path)

# Read in the metadata files and parse it for bulk samples of interest ---------
# Parsing code adapted from: https://github.com/AlexsLemonade/scpca-paper-figures/issues/98#issuecomment-2573164429
stopifnot(
  "Metadata files could not be found. Were they synced?" = all(file.exists(c(library_metadata_file, sample_metadata_file)))
)
library_metadata <- readr::read_tsv(library_metadata_file)
sample_metadata <- readr::read_tsv(sample_metadata_file)

# samples of interest
bulk_samples <- library_metadata |>
  dplyr::filter(seq_unit == "bulk") |>
  dplyr::pull(scpca_sample_id) |>
  unique()

# keep only solid tumors (N=142)
solid_bulk_samples_df <- sample_metadata |>
  dplyr::filter(scpca_sample_id %in% bulk_samples) |>
  dplyr::select(scpca_sample_id, diagnosis) |>
  dplyr::filter(
    !stringr::str_detect(diagnosis, "leukemia"),
    !stringr::str_detect(diagnosis, "Non-cancerous")
  )

# Create data frame to iterate over for syncing these bulk files
sync_bulk_df <- library_metadata |>
  dplyr::filter(
    scpca_sample_id %in% solid_bulk_samples_df$scpca_sample_id,
    seq_unit == "bulk"
  ) |>
  dplyr::mutate(
    output_dir = file.path(
      data_path,
      scpca_project_id, 
      scpca_sample_id,
      scpca_library_id
    )) |>
  # only keep columns needed for syncing
  dplyr::select(
    scpca_library_id, 
    output_dir
  )


# Sync each quant.sf file to its target directory -------
sync_bulk_df |>
  purrr::pwalk(
    \(scpca_library_id, output_dir) {
      quant_sf_s3 <- file.path(s3_path, scpca_library_id)
      
      # need trailing slash on s3 path
      sync_call <- glue::glue("aws s3 sync '{quant_sf_s3}/' '{output_dir}' --exclude '*' --include 'quant.sf'")
      system(sync_call)
  }
)
