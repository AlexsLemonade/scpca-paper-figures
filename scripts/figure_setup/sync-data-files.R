# This script is used to download any data needed to generate figures 
# All files will be saved in a folder named `s3_files`

renv::load()

# sync results files for SCPCP000001 
qc_s3_files <- "s3://nextflow-ccdl-results/scpca-prod/results/SCPCP000001/SCPCS000001"
qc_local_files <- here::here("s3_files", "SCPCS000001")
fs::dir_create(qc_local_files)

sync_call <- glue::glue("aws s3 cp '{qc_s3_files}' '{qc_local_files}' --exclude '*' --include '*.rds' --recursive")
system(sync_call)

# sync results for benchmarking libraries 
# define ids and s3 directories 
benchmarking_run_ids <- c("SCPCR000003", "SCPCR000126", "SCPCR000127", "SCPCR000220", "SCPCR000221", "SCPCR000495")
alevin_s3_dir <- "s3://nextflow-ccdl-results/scpca/alevin-fry-unfiltered-quant"
cellranger_s3_dir <- "s3://nextflow-ccdl-results/scpca/cellranger-quant"

# create local directory 
local_benchmark_dir <- here::here("s3_files", "benchmarking_results")
fs::dir_create(local_benchmark_dir)

# create includes statement to use for copying
aws_includes <- paste("--include '", benchmarking_run_ids, "*'", sep = '', collapse = ' ')

# sync alevin-fry 
sync_call <- glue::glue(
  "aws s3 cp '{alevin_s3_dir}' '{local_benchmark_dir}' --exclude '*' {aws_includes} --exclude '*.rad' --recursive"
)
system(sync_call)

# sync cell ranger 
# exclude large bam files 
sync_call <- glue::glue(
  "aws s3 cp '{cellranger_s3_dir}' '{local_benchmark_dir}' --exclude '*' {aws_includes} --exclude '*/SC_RNA_COUNTER_CS/*' --exclude '*/SPATIAL_RNA_COUNTER_CS/*' --exclude '*.bam' --exclude '*.bam.bai' --recursive"
)
system(sync_call)
