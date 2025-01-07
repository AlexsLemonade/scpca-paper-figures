# This script is used to download and save any reference files from S3 to a local directory
# All files will be saved in a folder named `s3_files/reference_files`

renv::load()

# create a directory to store S3 files ----------
local_s3_dir <- here::here("s3_files", "reference_files")
fs::dir_create(local_s3_dir)


# path on S3 to mito and t2g files -----
annotation_s3 <- "s3://scpca-references/homo_sapiens/ensembl-104/annotation"

# sync mito file -----
mito_s3 <- file.path(annotation_s3, "Homo_sapiens.GRCh38.104.mitogenes.txt")
sync_call <- glue::glue("aws s3 cp '{mito_s3}' '{local_s3_dir}'")
system(sync_call)

# sync t2g file -----
t2g_s3 <- file.path(annotation_s3, "Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv")
sync_call <- glue::glue("aws s3 cp '{t2g_s3}' '{local_s3_dir}'")
system(sync_call)


# sync SingleR models  ---------
singler_s3_dir <- "s3://scpca-references/celltype/singler_models"
singler_local_dir <- here::here("s3_files", "reference_files", "singler_models")
sync_call <- glue::glue("aws s3 cp '{singler_s3_dir}' '{singler_local_dir}' --recursive")
system(sync_call)
