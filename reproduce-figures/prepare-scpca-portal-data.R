#!/usr/bin/env Rscript

# This script organizes and parses data obtained from the ScPCA Portal so that figures and analyses can be reproduced.
# This script will create two directories you'll need:
# - `s3_files`: This directory will contain ScPCA data files needed to reproduce tables and figures.
# - `scpca-data`: This directory will contain ScPCA data files needed to reproduce the bulk expression analysis.
#
# For full details and instructions, please refer to the current directory's `README.md`.
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
# 3. In addition, you will need to download the merged version of project SCPCP000003, again in SingleCellExperiment format
# 4. Without renaming downloaded files, place all downloaded project ZIP files a single directory for input to this script.
#    Optionally, you can also store the merged version of SCPCP000003 in this directory, but it is not necessary.

# Then, you can run this script as follows:
# 
#   Rscript prepare-scpca-portal-data.R \
#     --portal_projects_dir <path to directory with all project zip files> \
#     --merged_sce_path <path to merged SCE file>
# 
# By default, the output directories will be saved in the location in this repository where code expects them:
# - `s3_files/` will be placed in the top level of the `scpca-paper-figures` repository
# - `scpca-data/` will be placed inside the bulk analysis data directory, `scpca-paper-figures/pseudobulk-bulk-prediction/data/`
# You can provide custom directories if you prefer; use the `--help` flag to see all script options.
# Note that this script will not overwrite existing output directories without the `--overwrite` flag.
#


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
    help = "Path to directory containing all project zip files downloaded from the ScPCA Portal"
  ),
  make_option(
    "--merged_sce_path",
    type = "character",
    help = "Path to the merged SCPCP000003 project file as a zip file"
  ),
  make_option(
    "--output_dir",
    type = "character",
    default = here::here(), 
    help = "Output directory where directories `s3_files` and `scpca-data` will be saved. By default, this is `./scpca-portal-data`."
  ),
  make_option(
    "--project_whitelist",
    type = "character",
    default = here::here("sample-info", "project-whitelist.txt"), 
    help = "Path to file containing all projects considered in the manuscript."
  ),  
  make_option(
    "--overwrite",
    action = "store_true",
    default = FALSE,
    help = "Whether to overwrite output files if they already exist in the provided output directory (default FALSE)"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))


# Genes whose expression should be kept in consensus cell type marker gene TSVs
marker_gene_ref_file <- "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/tags/v0.2.2/analyses/cell-type-consensus/references/validation-markers.tsv"

# Setup --------------------

# Check files and directories
stopifnot(
  "Portal projects directory not found" = dir.exists(opts$portal_projects_dir), 
  "Merged SCE for project SCPCP000003 not found" = file.exists(opts$merged_sce_path),
  "Output directory not specified" = !is.null(opts$output_dir), 
  "Project whitelist could not be found" = file.exists(opts$project_whitelist)
)

# Define and check output directories
output_dir_s3_files <- file.path(opts$output_dir, "s3_files")
output_dir_analysis_data <- file.path(opts$output_dir, "scpca_data")

if ((dir.exists(output_dir_s3_files) | dir.exists(output_dir_analysis_data)) & !opts$overwrite) {
  stop("Output directories already exist. To overwrite them, use the --overwrite flag.")
}

# Define additional nested directories and relevant sample/library ids
output_dir_consensus_files <- file.path(output_dir_s3_files, "cell-type-consensus-results")
output_dir_celltype_results <- file.path(output_dir_s3_files, "celltype_results")

# Create directories which will be assumed to exist later
fs::dir_create(output_dir_analysis_data)
fs::dir_create(output_dir_celltype_results)

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

# The samples' processed.rds files should be saved to output_dir_celltype_results
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

# These are the directories to save to output_dir_analysis_data
# They should contain processed files organized as expected _AND_ with the bulk quant TSV file
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


# Copy over the merged SCPCP000003 file ----------------
output_merged_sce_path <- file.path(output_dir_s3_files, "SCPCP000003", "SCPCP000003_merged.rds")
fs::dir_create(dirname(output_merged_sce_path))
fs::file_copy(
  opts$merged_sce_path, 
  output_merged_sce_path
)


# Map over project zip files to organize files for reproducibility ------------------
input_zips |>
  purrr::iwalk(
    \(project_zip, project_id){

      project_consensus_dir <- file.path(output_dir_consensus_files, project_id)
      
      # Unzip the directory, specifying the target in consensus
      unzip(project_zip, exdir = project_consensus_dir)

      # First, remove all the qc htmls
      system(glue::glue("rm {project_consensus_dir}/*/*html"))
      
      # Define all sample directories present in this project, _excluding_ multiplexed samples
      # This is because multiplexed samples aren't used in any figures created from ScPCA data files
      sample_dirs <- list.dirs(
        project_consensus_dir, 
        full.names = TRUE
      ) |>
        purrr::set_names(\(x) {basename(x)}) |>
        # get rid of the project folder & multiplexed samples
        purrr::discard_at(\(x) {stringr::str_detect(x, "SCPCP\\d{6}")}) |>
        purrr::discard_at(\(x) {stringr::str_detect(x, "_")}) 

      ####### Step 1: Copy RDS files to target directories #######
      sample_dirs |>
        purrr::iwalk(
          \(sample_dir, sample_id) {
            
            # Copy processed rds files to celltype_results_files
            # The next step will parse these files into TSVs
            if (sample_id %in% celltype_results_samples) {
              system(glue::glue("cp {sample_dir}/*processed.rds {output_dir_celltype_results}"))
            } 
            
            # Copy rds files to standalone sample dir
            if (sample_id %in% standalone_sample_dirs) {
              system(glue::glue("cp -r {sample_dir} {output_dir_s3_files}"))
            } 
          
            # We can now remove the (un)filtered files
            system(glue::glue("rm {sample_dir}/*filtered.rds"))
      })

      # If this project is used in the bulk analysis,copy project over as well
      if (project_id %in% bulk_projects) {
        fs::dir_copy(
          project_consensus_dir, 
          output_dir_analysis_data
        )
      }
      
      ####### Step 2: Prepare consensus cell type files #######
      sample_dirs |>
        purrr::walk(
          \(sample_dir) {
              
            # Define input and output files
            processed_rds_file <- list.files(
              path = sample_dir,
              pattern = "_processed\\.rds$", 
              full.names = TRUE
            )
            stopifnot("Expected a single processed rds file." = length(processed_rds_file) == 1)
            
            # Define output file names
            output_consensus_tsv <- stringr::str_replace(
              processed_rds_file, "_processed\\.rds$", "_processed_consensus-cell-types.tsv.gz"
            )
      
            output_markers_tsv <- stringr::str_replace(
              processed_rds_file, "_processed\\.rds$", "_processed_marker-gene-expression.tsv.gz"
            )    
            
            # Read in the processed SCE
            sce <- readRDS(processed_rds_file)
            
            # Prepare and export the consensus cell types table
            data.frame(
              project_id = metadata(sce)$project_id, 
              sample_id  = metadata(sce)$sample_id, 
              library_id = metadata(sce)$library_id, 
              sample_type = paste0(metadata(sce)$sample_type, collapse = ",")
              barcodes = colnames(sce), 
              # TODO TEMPORARY TO TEST THE CODE!!!!!!!!!!!
              consensus_annotation = sce$singler_celltype_annotation 
            ) |>
              # export
              readr::write_tsv(output_consensus_tsv)
            
            # Prepare and export the marker gene expression table
            # This code is adapted from the consensus cell type module in OpenScPCA-nf:
            # https://github.com/AlexsLemonade/OpenScPCA-nf/blob/a97e88c411cabb1b41e1fce6a7735311fb435051/modules/cell-type-consensus/resources/usr/bin/assign-consensus-celltypes.R#L205
            
            expressed_markers <- rowData(sce) |>
              as.data.frame() |>
              dplyr::filter(detected > 0) |>
              dplyr::pull(gene_ids) |>
              intersect(marker_genes)
      
            # get logcounts from sce for expressed genes, unless the length is 0 in which case fill NAs
            if (length(expressed_markers) == 0) {
              gene_exp_df <- data.frame(
                barcodes = rownames(sce),
                library_id = metadata(sce)$library_id,
                ensembl_gene_id = NA,
                logcounts = NA
              )
            } else {
              gene_exp_df <- scuttle::makePerCellDF(
                sce,
                features = expressed_markers,
                assay.type = "logcounts",
                use.coldata = "barcodes",
                use.dimred = FALSE
              ) |>
                tidyr::pivot_longer(starts_with("ENSG"), names_to = "ensembl_gene_id", values_to = "logcounts") |>
                # add library id column
                dplyr::mutate(library_id = metadata(sce)$library_id, .before = 0)
            }
            
            # export the marker gene expression tsv
            readr::write_tsv(gene_exp_df, output_markers_tsv)
            
            # Now that we have prepared consensus tsvs, we can remove the processed SCE itself
            fs::file_delete(processed_rds_file)
      })
})
