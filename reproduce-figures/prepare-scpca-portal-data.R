#!/usr/bin/env Rscript
#
# This script organizes and parses data obtained from the ScPCA Portal so that figures and analyses can be reproduced.
# This script will create three directories you'll need:
# - `s3_files`: This directory will contain ScPCA data, metadata, and reference files needed to reproduce figures, tables, and the bulk expression analysis.
# - `analysis/pseudobulk-bulk-analysis/data/scpca_data`: This directory will contain ScPCA data files needed to reproduce the bulk expression analysis.
# - `analysis/pseudobulk-bulk-analysis/data/references`: This directory will contain reference files needed to reproduce the bulk expression analysis.
#
#####################################################################################
### CAUTION! YOU WILL NEED AT LEAST 170 GB OF AVAILABLE SPACE TO RUN THIS SCRIPT. ###
###                                                                               ###
### In addition, it will take roughly 90 minutes to run this script.              ###
#####################################################################################
#
# For full details and usage instructions, including data you need to download before running this script,
# please refer to the current directory's `README.md`.
# Briefly, this script can be run as follows:
#
#
#   Rscript prepare-scpca-portal-data.R \
#     --portal_download_dir <path to directory with portal-wide download> \
#     --portal_metadata_path <path to portal-wide metadata zip file>

renv::load()
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# Parse options --------
option_list <- list(
  make_option(
    "--portal_download_dir",
    type = "character",
    help = "Path to directory containing the portal-wide download TODO: COMPRESSED OR NOT? WE'RE JUST GOING TO UNCOMPRESS IT ANYWAYS."
  ),
  make_option(
    "--portal_metadata_path",
    type = "character",
    help = "Path to the portal-wide metadata file as a compressed ZIP file"
  ),
  make_option(
    "--merged_project_path",
    type = "character",
    help = "Path to the the SCPCP000004 merged object as a compressed ZIP file"
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
    default = here::here("analysis", "pseudobulk-bulk-prediction", "data"),
    help = "Output directory for data used in the bulk RNA-Seq analysis. Default is `analysis/pseudobulk-bulk-prediction/data`, which is where code expects it to be."
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
marker_gene_ref_file <- "https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/refs/tags/v0.9.2/references/validation-markers.tsv"

# Functions used to prepare consensus cell types
utils_file <- here::here("reproduce-figures", "utils.R")
source(utils_file)

# Setup --------------------

# Check files and directories - for user-provided, first if they were specified, then if they exist
stopifnot(
  "Path to unzipped ScPCA Portla download must be specified with --portal_download_dir" = !is.null(opts$portal_download_dir),
  "Portal-wide metadata TSV must be specified with --portal_metadata_path" = !is.null(opts$portal_metadata_path)
)
stopifnot(
  "Portal projects directory not found" = dir.exists(opts$portal_download_dir),
  "Portal-wide metadata file not found" = file.exists(opts$portal_metadata_path),
  "Project whitelist could not be found" = file.exists(opts$project_whitelist)
)


# Check output directories

# first, set up the target bulk dirs
bulk_scpca_data_dir <- file.path(opts$bulk_data_dir, "scpca_data")
bulk_references_dir <- file.path(opts$bulk_data_dir, "references")
output_dirs <- c(
  opts$s3_files_dir,
  bulk_scpca_data_dir,
  bulk_references_dir
)
if (any(dir.exists(output_dirs))) {
  if (!opts$overwrite) {
    stop("Output directories already exist. To overwrite them, use the --overwrite flag.")
  } else {
    # use system to remove, in case any don't actually exist
    system(glue::glue("rm -rf {paste(output_dirs, collapse = ' ')}"))
  }
}

# Define and create directories
consensus_files_dir <- file.path(opts$s3_files_dir, "cell-type-consensus-results")
celltype_results_dir <- file.path(opts$s3_files_dir, "celltype_results")
s3_files_reference_dir <- file.path(opts$s3_files_dir, "reference_files")
SCPCP000003_sce_dir <- file.path(opts$s3_files_dir, "SCPCP000003")
SCPCP000004_sce_dir <- file.path(opts$s3_files_dir, "SCPCP000004")


# these directories are for temporary file stored in scratch
portal_metadata_scratch_dir <- file.path(opts$scratch_dir, "portal-metadata")
bulk_metadata_scratch_dir <- file.path(opts$scratch_dir, "bulk-metadata")
citeseq_metadata_scratch_dir <- file.path(opts$scratch_dir, "citeseq-project-metadata")


fs::dir_create(c(
  output_dirs,
  consensus_files_dir,
  celltype_results_dir,
  s3_files_reference_dir,
  SCPCP000003_sce_dir,
  SCPCP000004_sce_dir,
  bulk_metadata_scratch_dir,
  citeseq_metadata_scratch_dir,
  portal_metadata_scratch_dir
))


# Read input metadata file and define output metadata files
# Note that we need directory definitions before reading this in
# first, confirm it's a zip file
stopifnot(
  "The portal metadata file is expected to be a compressed ZIP file." = stringr::str_ends(tolower(opts$portal_metadata_path), "zip")
)
unzip(opts$portal_metadata_path, exdir = portal_metadata_scratch_dir)
portal_metadata <- readr::read_tsv(
  file.path(portal_metadata_scratch_dir, "metadata.tsv") # The downloaded metadata from the portal is literally named metadata.tsv
)
sample_metadata_file <- file.path(opts$s3_files_dir, "scpca-sample-metadata.tsv")
library_metadata_file <- file.path(opts$s3_files_dir, "scpca-library-metadata.tsv")


# These are all output directories at the top-level of s3_files
### Some need just processed, and some also need (un)filtered, but scripts specify which one
### We can therefore retain all RDS files when copying; extra files will not be used either way
standalone_sample_dirs <- c(
  "SCPCS000001",
  "SCPCS000216",
  "SCPCS000049",
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

# These are the directories to save to bulk_scpca_data_dir
# They should contain processed files organized as expected and the bulk quant TSV file
bulk_projects <- c(
  "SCPCP000001",
  "SCPCP000002",
  "SCPCP000006",
  "SCPCP000009",
  "SCPCP000017"
)

# Projects with CITE-seq technology which we'll assign technology to when preparing metadata
citeseq_projects <- c("SCPCP000003", "SCPCP000007", "SCPCP000008")

# Samples whose libraries we need to obtain for making the SCPCP000003 merged figure
SCPCP000003_samples <- c("SCPCS000050", "SCPCS000051", "SCPCS000053", "SCPCS000054")

# Sample to save to SCPCP000004, along with SCPCP000004_merged.rds
nb_sample <- "SCPCS000112"

# Define marker genes for recreating expression matrices used for consensus cell type dot plots
marker_genes <- readr::read_tsv(marker_gene_ref_file) |>
  dplyr::pull(ensembl_gene_id) |>
  unique()

# Define input files ------------------------

# Read whitelist so we can check that all projects are present
expected_projects <- readLines(opts$project_whitelist)

# Check for all project folders
expected_singlecell <- glue::glue("{expected_projects}_single-cell")
expected_bulk <- glue::glue("{bulk_projects}_bulk_rna")
stopifnot(
  "The provided portal download directory does not contain all expected single-cell projects. Please confirm the download was successful." = 
    all( dir.exists(file.path(opts$portal_download_dir, expected_singlecell)) ), 
  "The provided portal download directory does not contain all expected bulk projects. Please confirm the download was successful." = 
    all( dir.exists(file.path(opts$portal_download_dir, expected_bulk)) ) 
)

# Map over projects to organize files for reproducibility -------------
expected_projects |>
  purrr::walk(
    \(project_id) {
      
      project_download_dir <- file.path(opts$portal_download_dir, glue::glue("{project_id}_single-cell"))

      # Define all sample directories present in this project, _excluding_ multiplexed samples
      # This is because multiplexed samples aren't used in any figures created from ScPCA data files
      sample_dirs <- list.dirs(
        project_download_dir,
        full.names = TRUE
      ) |>
        purrr::set_names(\(x) {
          basename(x)
        }) |>
        # get rid of itself & multiplexed samples
        purrr::discard_at(\(x) {
          stringr::str_detect(x, project_id)
        }) |>
        purrr::discard_at(\(x) {
          stringr::str_detect(x, "_")
        })

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
          }
        )

      # If this project is SCPCP000003, copy over a few of its processed files
      if (project_id == "SCPCP000003") {
        SCPCP000003_samples |>
          purrr::map(
            \(sample_id) {
              # use system since we need to use a glob here
              system(
                glue::glue("cp {project_download_dir}/{sample_id}/*_processed.rds {SCPCP000003_sce_dir}")
              )
            }
          )
      }

      # If this project is SCPCP000004, copy over its sample and merged object
      if (project_id == "SCPCP000004") {
        fs::dir_copy(
          file.path(project_download_dir, nb_sample), 
          SCPCP000004_sce_dir, 
          overwrite = TRUE
        )
        
        fs::dir_copy(
          file.path(opts$merged_project_path), # TODO whats the zip situation?
          SCPCP000004_sce_dir, 
          overwrite = TRUE
        )
      }

  
      ####### Step 2: Copy bulk results as needed to target analysis directories #######
      # note that we can use these bulk TSVs for later metadata parsing
      if (project_id %in% bulk_projects) {
        fs::dir_copy(
          file.path(opts$portal_download_dir, glue::glue("{project_id}_bulk_rna")),
          file.path(file.path(bulk_scpca_data_dir, project_id)),
          overwrite = TRUE
        )
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

                # Read in the processed SCE
                sce <- readRDS(rds)

                # Prepare and export consensus tsv
                prepare_consensus_tsv(sce, output_consensus_tsv)

                # Prepare and export marker gene expression tsv
                prepare_gene_expression_tsv(sce, marker_genes, output_gene_expr_tsv)
              })
          }
        )
    }
  )


# Prepare sample and library metadata files ---------------------------

# Create and export sample metadata table
# NOTE: CONFRIMED WE STILL NEED THIS!!!!!!!!
# those 4 samples are only going to appear in the portal-wide metadata
prepare_sample_metadata(
  portal_metadata,
  sample_metadata_file
)

# Create and export library metadata table
prepare_library_metadata(
  portal_metadata,
  opts$portal_download_dir,
  library_metadata_file
)

# Now we can remove the remaining files in scratch
fs::dir_delete(c(
  bulk_metadata_scratch_dir,
  citeseq_metadata_scratch_dir
))

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
  glue::glue("aws s3 cp s3://scpca-references/celltype/singler_models {s3_files_reference_dir}/singler_models --recursive --exclude '*' --include '*.rds' --no-sign-request")
)

# Panglao marker gene references
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/cellassign_references/bone-and-soft-tissue_PanglaoDB_2020-03-27.tsv {bulk_references_dir} --no-sign-request")
)
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/cellassign_references/brain-compartment_PanglaoDB_2020-03-27.tsv {bulk_references_dir} --no-sign-request")
)
system(
  glue::glue("aws s3 cp s3://scpca-references/celltype/cellassign_references/kidney-compartment_PanglaoDB_2020-03-27.tsv {bulk_references_dir} --no-sign-request")
)
