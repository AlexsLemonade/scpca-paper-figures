# This script is used to download any data needed to generate figures 
# All files will be saved in a folder named `s3_files`

renv::load()

# sync results files for SCPCS000001 (for QC plots)
qc_s3_files <- "s3://nextflow-ccdl-results/scpca-prod/results/SCPCP000001/SCPCS000001"
qc_local_files <- here::here("s3_files", "SCPCS000001")
fs::dir_create(qc_local_files)

sync_call <- glue::glue("aws s3 cp '{qc_s3_files}' '{qc_local_files}' --exclude '*' --include '*.rds' --recursive")
system(sync_call)

# sync results for SCPCS000216 (for ADT plots)
adt_s3_files <- "s3://nextflow-ccdl-results/scpca-prod/results/SCPCP000007/SCPCS000216"
adt_local_files <- here::here("s3_files", "SCPCS000216")
fs::dir_create(adt_local_files)

sync_call <- glue::glue("op run -- aws s3 cp '{adt_s3_files}' '{adt_local_files}' --exclude '*' --include '*.rds' --recursive")
system(sync_call)
