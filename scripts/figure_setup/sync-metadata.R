# This script is used to download and save metadata from S3 to a local directory
# All files will be saved in a folder named `s3_files`

renv::load()

# path to metadata stored on S3
sample_info_s3 <- 's3://ccdl-scpca-data/sample_info'
sample_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-sample-metadata.tsv'
library_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
project_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-project-metadata.tsv'
celltype_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-project-celltype-metadata.tsv'

# project_whitelist
project_whitelist <- here::here("sample-info", "project-whitelist.txt")

# create a directory to store S3 files 
local_s3_dir <- here::here("s3_files")
fs::dir_create(local_s3_dir)

# Copy sample metadata 
local_sample_metadata <- file.path(local_s3_dir, "scpca-sample-metadata.tsv")
sync_call <- paste('aws s3 cp', sample_metadata_s3, local_s3_dir, sep = " ")
system(sync_call)

# copy library metadata
local_library_metadata <- file.path(local_s3_dir, "scpca-library-metadata.tsv")
sync_call <- paste('aws s3 cp', library_metadata_s3, local_s3_dir, sep = " ")
system(sync_call)


# copy project metadata
local_project_metadata <- file.path(local_s3_dir, "scpca-project-metadata.tsv")
sync_call <- paste('aws s3 cp', project_metadata_s3, local_s3_dir, sep = " ")
system(sync_call)

# grab project whitelist 
projects <- readLines(project_whitelist)

# read in project metadata, filter to whitelist projects, and grab metadata file names
project_metadata_files <- readr::read_tsv(local_project_metadata) |> 
  dplyr::filter(scpca_project_id %in% projects) |> 
  dplyr::pull(metadata_file)

# copy all individual project metaddata files 
local_project_dir <- file.path(local_s3_dir, "project-metadata")
fs::dir_create(local_project_dir)

# grab project metadata files from s3
glue::glue(
  "aws s3 sync '{sample_info_s3}' '{local_project_dir}' --exclude '*' --include '{project_metadata_files}' --exact-timestamps"
) |>
  purrr::walk(system)


# copy celltype ref metadata 
local_celltype_metadata <- file.path(local_s3_dir, "scpca-project-celltype-metadata.tsv")
sync_call <- glue::glue("aws s3 cp '{celltype_metadata_s3}' '{local_celltype_metadata}'")
system(sync_call)
