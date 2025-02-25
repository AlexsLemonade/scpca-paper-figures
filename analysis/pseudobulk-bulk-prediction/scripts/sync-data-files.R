# This script prepares files needed for analysis:
# 1. First, we identify samples of interest as non-multiplexed samples from solid tumors with paired bulk data
# 2. Second, we sync all associated bulk `quant.sf` files and single-cell `_processed.rds` files from S3 needed for analysis
# 3. Third, we sync all bulk counts files `_bulk_quant.tsv` from S3
# All files are organized in `<project id>/<sample id>/` (except bulk counts which are only per project)
# In addition, we sync PanglaoDB marker gene files used for cell type annotation for our projects of interest.
# Finally, we also export a TSV mapping gene symbols to ensembl ids, since ESTIMATE uses gene symbols


renv::load()
library(optparse)


# Parse options --------
option_list <- list(
  make_option(
    "--output_dir",
    type = "character",
    help = "Output directory to save synced ScPCA data files organized by project/sample"
  ),
  make_option(
    "--reference_dir",
    type = "character",
    help = "Output directory to save synced PanglaoDB cell type reference files"
  ),
  make_option(
    "--map_file",
    type = "character",
    help = "Path to output TSV file mapping bulk sample and library ids",
  ),
  make_option(
    "--sample_metadata_file",
    type = "character",
    default = here::here("s3_files", "scpca-sample-metadata.tsv"),
    help = "Path to ScPCA sample metadata file."
  ),
  make_option(
    "--library_metadata_file",
    type = "character",
    default = here::here("s3_files", "scpca-library-metadata.tsv"),
    help = "Path to ScPCA library metadata file."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Define files and directories ---------

stopifnot(
  "Metadata files could not be found. Were they synced?" = all(file.exists(c(opts$library_metadata_file, opts$sample_metadata_file)))
)
s3_bulk_tpm_path <- "s3://nextflow-ccdl-results/scpca-prod/checkpoints/salmon" #/library_id/quant.sf
s3_processed_path <- "s3://nextflow-ccdl-results/scpca-prod/results" #/project_id/sample_id/library_id_processed.rds AND project_id/bulk/<project_id>_bulk_quant.tsv

fs::dir_create(opts$output_dir)
fs::dir_create(opts$reference_dir)


## Determine samples of interest -----------

library_metadata <- readr::read_tsv(opts$library_metadata_file, show_col_types = FALSE)
sample_metadata <- readr::read_tsv(opts$sample_metadata_file, show_col_types = FALSE)

# we can't include multiplexed in this analysis, so we'll find those for removal
cellhash_samples <- library_metadata |>
  dplyr::filter(stringr::str_detect(technology, "cellhash")) |>
  dplyr::pull(scpca_sample_id) |>
  stringr::str_split(pattern = ";") |>
  purrr::reduce(union)

# all possible samples of interest
all_bulk_samples <- library_metadata |>
  dplyr::filter(
    !(scpca_sample_id %in% cellhash_samples),
    seq_unit == "bulk"
  ) |>
  dplyr::pull(scpca_sample_id) |>
  unique()


# keep only solid tumors with directly paired single-cell and bulk
solid_bulk_samples <- sample_metadata |>
  dplyr::filter(scpca_sample_id %in% all_bulk_samples) |>
  dplyr::filter(
    !stringr::str_detect(diagnosis, "leukemia"),
    !stringr::str_detect(diagnosis, "Non-cancerous")
  ) |>
  dplyr::pull(scpca_sample_id)

# Create data frame to iterate over for syncing single-cell files
# we do this first since there are fewer single-cell than corresponding bulk
sync_sc_df <- library_metadata |>
  dplyr::filter(
    scpca_sample_id %in% solid_bulk_samples,
    seq_unit %in% c("cell", "nucleus")
  ) |>
  dplyr::mutate(
    s3_dir = file.path(s3_processed_path, scpca_project_id, scpca_sample_id),
    output_dir = file.path(opts$output_dir, scpca_project_id, scpca_sample_id)
  ) |>
  # only keep columns needed for syncing or checking
  dplyr::select(
    scpca_sample_id,
    scpca_library_id,
    output_dir,
    s3_dir
  )


# Create data frame to iterate over for syncing bulk TPM files
sync_bulk_tpm_df <- library_metadata |>
  dplyr::filter(
    scpca_sample_id %in% sync_sc_df$scpca_sample_id,
    seq_unit == "bulk"
  ) |>
  dplyr::mutate(
    s3_dir = file.path(s3_bulk_tpm_path, scpca_library_id),
    output_dir = file.path(opts$output_dir, scpca_project_id, scpca_sample_id)
  ) |>
  # only keep columns needed for syncing or checking
  dplyr::select(
    scpca_sample_id, # this column is used for checking, not syncing
    output_dir,
    s3_dir
  )



# Create data frame to iterate over for syncing bulk counts TSV files
sync_bulk_counts_df <- library_metadata |>
  dplyr::filter(
    scpca_sample_id %in% sync_sc_df$scpca_sample_id,
    seq_unit == "bulk"
  ) |>
  dplyr::mutate(
    s3_dir = file.path(s3_processed_path, scpca_project_id, "bulk"),
    output_dir = file.path(opts$output_dir, scpca_project_id)
  ) |>
  # only keep columns needed for syncing
  dplyr::select(
    scpca_project_id,
    output_dir,
    s3_dir
  ) |>
  dplyr::distinct()


# Confirm we're matching all around
stopifnot(
  "Bulk and single-cell samples don't match" = setequal(sync_sc_df$scpca_sample_id, sync_bulk_tpm_df$scpca_sample_id)
)



# Sync quant.sf files ------------------
sync_bulk_tpm_df |>
  purrr::pwalk(
    \(scpca_sample_id, output_dir, s3_dir) {
      # need trailing slash on s3 path
      sync_call <- glue::glue("aws s3 sync '{s3_dir}/' '{output_dir}' --exclude '*' --include 'quant.sf'")
      system(sync_call)
    }
  )


# Sync the _processed.rds files ------------------
sync_sc_df |>
  purrr::pwalk(
    \(scpca_sample_id, scpca_library_id, output_dir, s3_dir) {
      sync_call <- glue::glue("aws s3 sync '{s3_dir}/' '{output_dir}' --exclude '*' --include '{scpca_library_id}_processed.rds'")
      system(sync_call)
  }
)

# Sync the bulk counts TSV files ------------------
sync_bulk_counts_df |>
  purrr::pwalk(
    \(scpca_project_id, output_dir, s3_dir) {
      sync_call <- glue::glue("aws s3 sync '{s3_dir}/' '{output_dir}' --exclude '*' --include '{scpca_project_id}_bulk_quant.tsv'")
      system(sync_call)
  }
)


# Sync Panglao gene sets -------------------------
# from https://github.com/AlexsLemonade/ScPCA-admin/blob/a2b2daf0caac336f391bac974bc1b808f29d5853/sample_info/scpca-project-celltype-metadata.tsv, we need:
# - `brain-compartment` for SCPCP000001, SCPCP000002, and SCPCP000009
# - `kidney-compartment` for SCPCP000006
# - `bone-and-soft-tissue` for SCPCP000017

c("brain-compartment", "kidney-compartment", "bone-and-soft-tissue") |>
  purrr::walk(
    \(compartment) {
      sync_call <- glue::glue("aws s3 sync 's3://scpca-references/celltype/cellassign_references/' '{opts$reference_dir}' --exclude '*' --include '{compartment}_PanglaoDB_2020-03-27.tsv'")
      system(sync_call)
    }
  )


# Finally, export a TSV of matching sample and library ids for bulk ------------
# we do this for later mapping because bulk counts TSV are labeled by library, not sample, ids
bulk_libraries_samples <- library_metadata |>
  dplyr::filter(
    scpca_sample_id %in% sync_sc_df$scpca_sample_id,
    seq_unit == "bulk"
  ) |>
  dplyr::select(scpca_sample_id, scpca_library_id)
readr::write_tsv(bulk_libraries_samples, opts$map_file)
