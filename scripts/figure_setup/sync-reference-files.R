# This script is used to download and save any reference files from S3 to a local directory
# All files will be saved in a folder named `s3_files/reference_files`

renv::load()

# path to mito file 
mito_s3 <- 's3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt'

# create a directory to store S3 files 
local_s3_dir <- here::here("s3_files", "reference_files")
fs::dir_create(local_s3_dir)

# Copy mito file
local_mito_file <- file.path(local_s3_dir, "Homo_sapiens.GRCh38.104.mitogenes.txt")
sync_call <- glue::glue("aws s3 cp '{mito_s3}' '{local_mito_file}'")
system(sync_call)
