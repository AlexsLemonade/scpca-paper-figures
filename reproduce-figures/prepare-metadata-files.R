# This script can be used to create versions of `scpca-sample-metadata.tsv` and `scpca-library-metadata.tsv` needed to
# reproduce manuscript figures, tables, and analyses.
#
# To run this script, you will first need to download all sample data from the ScPCA Portal: <https://scpca.alexslemonade.org/>
# From this page, click the `Get All Sample Metadata` button on the top-right of the page.
# Provide this TSV file as input to this script as:
#
# Rscript prepare-metadata-files.R --portal_sample_metadata <path to metadata.tsv>
#
# Files will be exported to a directory called `metadata-files/` by default but this can be specified with `--output_dir`.
# If metadata files are already present in the output directory, this script will not overwrite them by default.
# To overwrite existing metadata files, use the `--overwrite` flag.


renv::load()
library(optparse)

# Parse options --------
option_list <- list(
  make_option(
    "--portal_sample_metadata",
    type = "character",
    help = "Path to TSV of all sample metadata downloaded from ScPCA Portal"
  ),
  make_option(
    "--output_dir",
    type = "character",
    default = "metadata-files",
    help = "Path to output directory to save TSV files to. Default is a directory `metadata-files/`.
    If metadata files already exist in the given, this script will not overwrite them unless you specify the `--overwrite` flag."
  ),
  make_option(
    "--overwrite",
    action = "store_true",
    default = FALSE,
    help = "Whether to overwrite existing metadata files found in the output directory."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check that the input file exists
stopifnot("Portal sample metadata TSV file not found" = file.exists(opts$portal_sample_metadata))

# Create output directory
fs::dir_create(opts$output_dir)

# Define output file paths
sample_metadata_tsv <- file.path(opts$output_dir, "scpca-sample-metadata.tsv")
library_metadata_tsv <- file.path(opts$output_dir, "scpca-library-metadata.tsv")

# If files already exist, stop unless `--overwrite` flag was provided
if ((file.exists(sample_metadata_tsv) | file.exists(library_metadata_tsv)) & (!opts$overwrite)) {
  print("Metadata files already exist in the output directory. To overwrite them, use the --overwrite flag.")
  stop()
}

# Read input metadata file
portal_metadata <- readr::read_tsv(opts$portal_sample_metadata)

# Create sample metadata table
sample_metadata <- portal_metadata |>
  dplyr::select(scpca_project_id, scpca_sample_id, diagnosis, disease_timing, is_cell_line) |>
  # remove duplicate rows, which occur when there are multiple libraries per sample
  dplyr::distinct()

# Create library metadata table
library_metadata <- portal_metadata |>
  dplyr::select(scpca_project_id, scpca_sample_id, scpca_library_id, seq_unit, technology)

# Export tsv files
readr::write_tsv(sample_metadata, sample_metadata_tsv)
readr::write_tsv(library_metadata, library_metadata_tsv)
