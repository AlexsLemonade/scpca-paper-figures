#!/usr/bin/env Rscript

# This script organizes and parses data obtained from the ScPCA Portal so that figures and analyses can be reproduced.
# This script will create two directories you'll need:
# - `s3_files`: This directory will contain ScPCA data files needed to reproduce tables and figures.
# - `scpca-data`: This directory will contain ScPCA data files needed to reproduce the bulk expression analysis.
#
# For full details and instructions, please refer to the current directory's `README.md`.
#
#####################################################################################
### CAUTION! YOU WILL NEED AT LEAST 170 GB OF AVAILABLE SPACE TO RUN THIS SCRIPT. ###
#####################################################################################
#
#
# ########################### Instructions ############################
#
# Before running this script, you will need to take the following steps:
#
# 1. Ensure you have set up the `renv` environment in this project.
#    For more information, refer to the top-level `README.md` file.
# 2. Download all projects in the ScPCA Portal that are listed in the file `../sample-info/project-whitelist.txt`, with the following download options:
#   - Ensure the selected modality is Single-cell
#   - Ensure the selected Data Format is SingleCellExperiment (r)
#   - Do _not_ click to merge all objects into the same sample
#   - If the project contains multiplexed samples, you can exclude those samples (but it will not affect the code if they are included)
# 3. Download the merged version of project SCPCP000003 from the ScPCA Portal, again in SingleCellExperiment format
# 4. Download the portal-wide metadata from the ScPCA Portal using the "Get All Sample Metadata" button on the top-right of the portal homepage
# 5. Place all downlaoded project ZIP files (step 2) a single directory for input to this script.
#    Optionally, you can also store the merged version of SCPCP000003 and the portal metadata TSV in this directory, but it is not necessary.

# Then, you can run this script as follows:
#
#   Rscript prepare-scpca-portal-data.R \
#     --portal_projects_dir <path to directory with all project zip files> \
#     --merged_sce_path <path to merged SCE file> \
#     --portal_metadata_path <path to portal-wide metadata TSV>
#
# By default, the output directories will be saved in the location in this repository where code expects them:
# - `s3_files/` will be placed in the top level of the `scpca-paper-figures` repository
# - `scpca-data/` will be placed inside the bulk analysis data directory, `scpca-paper-figures/pseudobulk-bulk-prediction/data/`
# You can provide custom directories if you prefer; use the `--help` flag to see all script options.
# Note that this script will not overwrite existing output directories without the `--overwrite` flag.



renv::load()
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})


# Parse options --------
option_list <- list(
  make_option(
    "--portal_projects_dir",
    type = "character",
    help = "Path to directory containing project ZIP files downloaded from the ScPCA Portal"
  ),
  make_option(
    "--portal_metadata_path",
    type = "character",
    help = "Path to the the portal-wide metadata TSV"
  ),
  make_option(
    "--merged_sce_path",
    type = "character",
    help = "Path to the merged SCPCP000003 project ZIP file"
  ),
  make_option(
    "--s3_files_dir",
    type = "character",
    default = here::here("s3_files"),
    help = "Output directory for `s3_files`. Default is in the top-level of the repository, which is where code expects it to be."
  ),
  make_option(
    "--bulk_data_dir",
    type = "character",
    default = here::here("analysis", "pseudobulk-bulk-prediction", "data", "scpca-data"),
    help = "Output directory for data used in the bulk RNA-Seq analysis. Default is `analysis/pseudobulk-bulk-prediction/data/scpca-data`, which is where code expects it to be."
  ),
  make_option(
     "--scratch_dir",
    type = "character",
    default = here::here("reproduce-figures", "scratch"),
    help = "Scratch directory for holding temporary files during processing"
  ),
  make_option(
    "--overwrite",
    action = "store_true",
    default = FALSE,
    help = "Whether to overwrite output files if they already exist in the provided output directory (default FALSE)"
  ),
  make_option(
    "--project_whitelist",
    type = "character",
    default = here::here("sample-info", "project-whitelist.txt"),
    help = "Path to file containing all projects considered in the manuscript."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Genes whose expression should be kept in consensus cell type marker gene TSVs
marker_gene_ref_file <- "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/tags/v0.2.2/analyses/cell-type-consensus/references/validation-markers.tsv"

# Functions used to prepare consensus cell types
utils_file <- here::here("reproduce-figures", "utils.R")
source(utils_file)

# Setup --------------------

# Check files and directories - for user-provided, first if they were specified, then if they exist
stopifnot(
  "Path to directory with project ZIP files be specified with --portal_projects_dir" = !is.null(opts$portal_projects_dir),
  "Path to merged SCE for project SCPCP000003 must be specified with --merged_sce_path" = !is.null(opts$merged_sce_path),
  "Portal-wide metadata TSV must be specified with --portal_metadata_path" = !is.null(opts$portal_metadata_path)
)
stopifnot(
  "Portal projects directory not found" = dir.exists(opts$portal_projects_dir),
  "Merged SCE for project SCPCP000003 not found" = file.exists(opts$merged_sce_path),
  "Portal-wide metadata TSV not found" = file.exists(opts$portal_metadata_path),
  "Project whitelist could not be found" = file.exists(opts$project_whitelist)
)

# Check output directories
if ((dir.exists(opts$s3_files_dir) | dir.exists(opts$bulk_data_dir))) {
  if (!opts$overwrite) {
    stop("Output directories already exist. To overwrite them, use the --overwrite flag.")
  } else {
    # use system to remove, in case one doesn't exist since we used an or above
    system(glue::glue("rm -rf {opts$s3_files_dir} {opts$bulk_data_dir}"))
  }
}

# Define and create directories
consensus_files_dir   <- file.path(opts$s3_files_dir, "cell-type-consensus-results")
celltype_results_dir  <- file.path(opts$s3_files_dir, "celltype_results")
s3_files_reference_dir <- file.path(opts$s3_files_dir, "reference_files")
bulk_reference_dir    <- file.path(opts$bulk_data_dir, "references")

fs::dir_create(c(
  opts$scratch_dir,
  opts$bulk_data_dir,
  consensus_files_dir,
  celltype_results_dir,
  s3_files_reference_dir,
  bulk_reference_dir
))

# Define marker genes for recreating expression matrices used for consensus cell type dot plots
marker_genes <- readr::read_tsv(marker_gene_ref_file) |>
  dplyr::pull(ensembl_gene_id) |>
  unique()

# These are all output directories at the top-level of s3_files
### Some need just processed, and some also need (un)filtered, but scripts specify which one
### We can therefore retain all RDS files when copying; extra files will not be used either way
standalone_sample_dirs <- c(
  "SCPCS000001",
  "SCPCS000216",
  "SCPCS000251",
  "SCPCS000264"
)

# These samples' processed.rds files should be saved to celltype_results_dir
celltype_results_samples <- c(
  "SCPCS000001",
  "SCPCS000002",
  "SCPCS000004",
  "SCPCS000252",
  "SCPCS000258",
  "SCPCS000262",
  "SCPCS000222",
  "SCPCS000223",
  "SCPCS000224"
)

# These are the directories to save to opts$bulk_files_dir
# They should contain processed files organized as expected and the bulk quant TSV file
bulk_projects <- c(
  "SCPCP000001",
  "SCPCP000002",
  "SCPCP000006",
  "SCPCP000009",
  "SCPCP000017"
)


# Define input files ------------------------

# Read whitelist so we can check that all projects are present
expected_projects <- readLines(opts$project_whitelist)

# Define and check input projects
input_zips <- list.files(
  opts$portal_projects_dir,
  # allow for a token after .zip
  pattern = "^SCPCP\\d{6}_single-cell.+\\.zip\\?*.*",
  full.names = TRUE,
  ignore.case = TRUE # case insensitive regex
) |>
  # remove any merged objects
  purrr::discard(\(x){grepl("_merged_", x, ignore.case = TRUE)}) |>
  # name by project
  purrr::set_names(
    \(x) {stringr::str_split_i(basename(x), "_", 1)}
  )
stopifnot(
  "The provided input directory does not contain all expected projects. Zip files for all projects listed in the `project-whitelist.txt` file should be present." =
    setequal(names(input_zips), expected_projects)
)

# First, unzip and copy over the merged SCPCP000003 file -----------------------
merged_sce_dir <- file.path(opts$s3_files_dir, "SCPCP000003")
merge_scratch_dir <- file.path(opts$scratch_dir, "SCPCP000003_merged")
fs::dir_create(c(merged_sce_dir, merge_scratch_dir))

unzip(opts$merged_sce_path, exdir = merge_scratch_dir)
fs::file_move(
  file.path(merge_scratch_dir, "SCPCP000003_merged.rds"),
  merged_sce_dir
)

# clean up
fs::dir_delete(merge_scratch_dir)

# Map over project zip files to organize files for reproducibility -------------
input_zips |>
  purrr::iwalk(
    \(project_zip, project_id){

      # scratch directory to store unzipped files during processing
      project_scratch_dir <- file.path(opts$scratch_dir, project_id)
      fs::dir_create(project_scratch_dir)

      # Unzip the directory to scratch
      unzip(project_zip, exdir = project_scratch_dir)

      # Define all sample directories present in this project, _excluding_ multiplexed samples
      # This is because multiplexed samples aren't used in any figures created from ScPCA data files
      sample_dirs <- list.dirs(
        project_scratch_dir,
        full.names = TRUE
      ) |>
        purrr::set_names(\(x) {basename(x)}) |>
        # get rid of itself & multiplexed samples
        purrr::discard_at(\(x) {stringr::str_detect(x, project_id)}) |>
        purrr::discard_at(\(x) {stringr::str_detect(x, "_")})

      ####### Step 1: Copy RDS files to target directories #######
      sample_dirs |>
        purrr::iwalk(
          \(sample_dir, sample_id) {

            # Copy processed rds files to celltype_results_files
            # The next step will parse these files into TSVs
            if (sample_id %in% celltype_results_samples) {
              system(glue::glue("cp {sample_dir}/*processed.rds {celltype_results_dir}"))
            }

            # Copy rds files to standalone sample dir
            if (sample_id %in% standalone_sample_dirs) {
              standalone_dir <- file.path(opts$s3_files_dir, sample_id)
              fs::dir_create(standalone_dir)
              system(glue::glue("cp -r {sample_dir}/*rds {standalone_dir}"))
            }

      })

      # If this project is used in the bulk analysis, copy to bulk_project_dir
      if (project_id %in% bulk_projects) {
        # First, we'll remove the files that aren't needed to save disk space
        system(glue::glue("rm {project_scratch_dir}/*/*filtered.rds"))
        system(glue::glue("rm {project_scratch_dir}/*/*html"))

        fs::dir_copy(project_scratch_dir, opts$bulk_data_dir)
      }

      ####### Step 2: Prepare consensus cell type files #######
      sample_dirs |>
        purrr::iwalk(
          \(sample_dir, sample_id) {

            # Map over all rds files for each sample
            list.files(
              path = sample_dir,
              pattern = "_processed\\.rds$",
              full.names = TRUE
            ) |>
              purrr::walk(\(rds) {

                # Define output file paths
                outdir <- file.path(consensus_files_dir, project_id, sample_id)
                fs::dir_create(outdir)

                output_consensus_tsv <- file.path(
                  outdir,
                  stringr::str_replace(basename(rds), "_processed\\.rds$", "_processed_consensus-cell-types.tsv.gz")
                )
                output_gene_expr_tsv <- file.path(
                  outdir,
                  stringr::str_replace(basename(rds), "_processed\\.rds$", "_processed_marker-gene-expression.tsv.gz")
                )

                # Read in the SCE
                sce <- readRDS(rds)

                # Prepare and export consensus tsv
                prepare_consensus_tsv(sce, output_consensus_tsv)

                # Prepare and export marker gene expression tsv
                prepare_gene_expression_tsv(sce, marker_genes, output_gene_expr_tsv)

            })
      })

      # Now that this project has been processed, we can clean out the scratch directory
      fs::dir_delete(project_scratch_dir)
})


# Prepare sample and library metadata files ---------------------------

# Read input metadata file
portal_metadata <- readr::read_tsv(opts$portal_metadata_path)
sample_metadata_file <- file.path(opts$s3_files_dir, "scpca-sample-metadata.tsv")
library_metadata_file <- file.path(opts$s3_files_dir, "library-sample-metadata.tsv")

# Create and export sample metadata table
 portal_metadata |>
  dplyr::select(scpca_project_id, scpca_sample_id, diagnosis, disease_timing, is_cell_line) |>
  # remove duplicate rows, which occur when there are multiple libraries per sample
  dplyr::distinct() |>
  readr::write_tsv(sample_metadata_file)

# Create and export library metadata table
library_metadata <- portal_metadata |>
  dplyr::select(scpca_project_id, scpca_sample_id, scpca_library_id, seq_unit, technology) |>
  readr::write_tsv(library_metadata_file)


# Sync reference files from S3 ------------------------------------

# annotation files
system(
  glue::glue("aws s3 cp s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv {s3_files_reference_dir} --no-sign-request")
)
system(
  glue::glue("aws s3 cp s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt {s3_files_reference_dir} --no-sign-request")
)

# SingleR model files
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/singler_models/ {s3_files_reference_dir} --recursive --exclude '*' --include '*.rds' --no-sign-request")
)

# Panglao marker gene references
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/cellassign_references/bone-and-soft-tissue_PanglaoDB_2020-03-27.tsv {bulk_reference_dir} --no-sign-request")
)
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/cellassign_references/brain-compartment_PanglaoDB_2020-03-27.tsv {bulk_reference_dir} --no-sign-request")
)
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/cellassign_references/kidney-compartment_PanglaoDB_2020-03-27.tsv {bulk_reference_dir} --no-sign-request")
)
