
#' Get whitelist of samples to include in plots 
#'
#' @param project_metadata_file Path to file containing ScPCA project metadata 
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