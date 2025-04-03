#!/usr/bin/env Rscript


renv::load()
library(optparse)

# Parse options --------
option_list <- list(
  make_option(
    "--portal_projects_dir",
    type = "character",
    default = "/Users/sjspielman/Downloads/projects",
    help = "Path to directory containing all projects downloaded from the ScPCA Portal"
  ),
  make_option(
    "--merged_sce_path",
    type = "character",
    default = "/Users/sjspielman/Downloads/SCPCP000003_SINGLE-CELL_SINGLE-CELL-EXPERIMENT_MERGED_2025-04-03.zip",
    help = "Path to Portal downloaded file merged SCPCP000003 as a zip file"
  ),
  make_option(
    "--output_dir",
    type = "character",
    default = "/Users/sjspielman/Desktop/parsed-portal",
    help = "Output directory where directories `s3_files` and `scpca_data` will be saved"
  ),
  make_option(
    "--overwrite",
    action = "store_true",
    default = FALSE,
    help = "Whether to overwrite existing directories `s3_files` and `scpca_data` if they exist in the provided output directory (default FALSE)"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))


# Setup --------------------

# Check files and directories
stopifnot(
  "Portal projects directory not found" = dir.exists(opts$portal_projects_dir), 
  "Merged SCE for project SCPCP000003 not found" = file.exists(opts$merged_sce_path),
  "Output directory not specified" = !is.null(opts$output_dir)
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
bulk_projects <- c("SCPCP000001", "SCPCP000002", "SCPCP000006", "SCPCP000009", "SCPCP000017")

# Number of projects expected
n_projects <- 23

# Prepare files ------------------------

# Define and check input projects
input_zips <- list.files(
  opts$portal_projects_dir, 
  pattern = "^SCPCP\\d{6}_SINGLE-CELL_SINGLE-CELL-EXPERIMENT_\\d+{4}-\\d+{2}-\\d+{2}\\.zip",
  full.names = TRUE
) |>
  purrr::set_names(
    \(x) {stringr::str_split_i(basename(x), "_", 1)}
  )
stopifnot(
  "Not all projects were found in the provided input directory. There should be 23 projects total." = 
    length(input_zips) == n_projects
)


input_zips |>
  purrr::iwalk(
    \(project_zip, project_id){
      
      # temporary variables used for testing
      #project_id <- "SCPCP000001"
      #project_zip <- "/Users/sjspielman/Downloads/projects/SCPCP000001_SINGLE-CELL_SINGLE-CELL-EXPERIMENT_2025-04-03.zip"
      
      project_consensus_dir <- file.path(output_dir_consensus_files, project_id)
      
      # Unzip the directory, specifying the target in consensus
      unzip(project_zip, exdir = project_consensus_dir)

      # Define all sample directories present in this project
      sample_dirs <- list.dirs(
        project_consensus_dir, 
        full.names = TRUE
      ) |>
        purrr::set_names(
          \(x) {basename(x)}
        )
      sample_dirs <- sample_dirs[stringr::str_detect(names(sample_dirs), "SCPCS\\d{6}")]
      
      # First, remove all the qc htmls
      system(glue::glue("rm {project_consensus_dir}/*/*html"))

      ####### Step 1: Copy RDS files to target directories #######
      sample_dirs |>
        purrr::iwalk(
          \(sample_dir, sample_id) {
            
            # First, copy RDS files to target destinations as needed
            
            # Copy processed rds files to celltype_results_files
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

      # If this is a bulk project, copy project over as well
      if (project_id %in% bulk_projects) {
        fs::dir_copy(
          project_consensus_dir, 
          output_dir_analysis_data
        )
      }
      
      ####### Step 2: Prepare consensus cell type files #######
      #sample_dirs |>
      #  purrr::iwalk(
      #    \(sample_dir, sample_id) {
      #      # forthcoming!
      #  })
      
})
