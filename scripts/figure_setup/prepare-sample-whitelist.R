# This script establishes a sample whitelist for determining which 
#  samples to include in manuscript figures.

renv::load()

# Functions -----------

#' Get whitelist of samples to include in plots 
#'
#' @param project_metadata_file Path to file containing ScPCA project metadata 
#' @param project_whitelist A vector of project ID's 
#'
#' @return List of sample IDs that are present on the portal 
get_sample_whitelist <- function(project_metadata_file, project_whitelist) {
  
  # read in project metadata files and create sample whitelist 
  project_metadata_files <- readr::read_tsv(project_metadata_file) |> 
    dplyr::filter(scpca_project_id %in% project_whitelist) |> 
    dplyr::pull(metadata_file)
  # get full file paths to each project metadata file
  project_metadata_files <- here::here("s3_files", "project-metadata", project_metadata_files)
  
  # grab samples that are on the portal and create a whitelist 
  sample_whitelist <- project_metadata_files |> 
    purrr::map(\(file){
      sample_list <- readr::read_tsv(file) |> 
        dplyr::filter(on_portal) |> 
        dplyr::pull(scpca_sample_id)
      
      return(sample_list)
    }) |> 
    unlist() 
  
  return(sample_whitelist)
}


# Define file paths ---------------
project_metadata_file <- here::here("s3_files", "scpca-project-metadata.tsv")
sample_whitelist_file <- here::here("sample-info", "sample-whitelist.txt")
project_whitelist_file <- here::here("sample-info", "project-whitelist.txt")


# Create and export sample whitelist --------------
project_whitelist <- readLines(project_whitelist_file)

sample_whitelist <- get_sample_whitelist(project_metadata_file, project_whitelist)
readr::write_lines(sample_whitelist, sample_whitelist_file)

