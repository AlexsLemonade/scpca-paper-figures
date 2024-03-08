# This script is used to download and save metadata from S3 to a local directory
# All files will be saved in a folder named `s3_files`

# path to metadata stored on S3
sample_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-sample-metadata.tsv'
library_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
celltype_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-project-celltype-metadata.tsv'

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

# copy celltype ref metadata 
local_celltype_metadata <- file.path(local_s3_dir, "scpca-project-celltype-metadata.tsv")
sync_call <- glue::glue("aws s3 cp '{celltype_metadata_s3}' '{local_celltype_metadata}'")
system(sync_call)
