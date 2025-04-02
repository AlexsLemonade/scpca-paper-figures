#!/usr/bin/env Rscript

# This script is used to prepare consensus cell type file result files for upload to Zenodo
# Using files in `s3_files/cell-type-consensus-results`, this script subsets files to only columns used in the code, 
#  as well as columns needed to identify samples and exports files in the same directory structure.
# This script is only intended for use by the Data Lab.
# Before running this script, you will to run `sync-consensus-celltype-results.R`, using `op run --` as needed to manage credentials.
#
# Usage:
#
# Rscript prepare-consensus-cell-types-zenodo.R --input_dir <path to consensus cell type results> --output_dir <path to new directory with results for Zenodo upload>
#
# By default, the output directory is `./cell-type-consensus-results`. The output directory will be compressed as `.tar.gz`.

library(optparse)


# Parse options --------
option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    default = here::here("s3_files", "cell-type-consensus-results"),
    help = "Path to full consensus cell type results from OpenScPCA"
  ),
  make_option(
    "--output_dir",
    type = "character",
    default = "./cell-type-consensus-results",
    help = "Path to output directory to save consensus cell types intended for Zenodo upload"
  ), 
  make_option(
    "--overwrite",
    action = "store_true",
    default = FALSE,
    help = "If the output directory already exists, this flag controls whether to overwrite it"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# check input directoru
stopifnot("Input directory not found" = dir.exists(opts$input_dir))

# if --overwrite was specified, remove output directory if it exists
if (dir.exists(opts$output_dir)) {
  if (opts$overwrite) {
    fs::dir_delete(opts$output_dir)
  } else {
    print("Output directory already exists. Use --overwrite to overwrite it.")
    stop()
  } 
}

# create output directory
fs::dir_create(opts$output_dir)

# Define input paths to consensus cell type files
celltype_files <- list.files(
  opts$input_dir,
  pattern = "_processed_consensus-cell-types\\.tsv\\.gz$",
  full.names = TRUE,
  recursive = TRUE
)

# Define vector of columns to retain
retain_columns <- c(
  "project_id", 
  "sample_id", 
  "library_id", 
  "barcodes",
  "sample_type", 
  "consensus_annotation"
)

# Export subsetted tsvs with only columns of interest
celltype_files |> 
  purrr::walk(
    \(celltype_file) {
      
      # define output file
      output_file <- stringr::str_replace(
         celltype_file, 
         opts$input_dir, 
         opts$output_dir
      )
      
      # create output directory
      fs::dir_create(dirname(output_file))

      # subset and export data frame
      readr::read_tsv(celltype_file) |>
        dplyr::select(all_of(retain_columns)) |>
        readr::write_tsv(output_file)
    }
)

# Finally, compress the output directory and remove the uncompressed directory
tar(
  tarfile = glue::glue(
    "{stringr::str_remove(opts$output_dir, '%')}.tar.gz"
  ), 
  files = opts$output_dir, 
  compression = "gz"
)
fs::dir_delete(opts$output_dir)
