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
af_local_dir <- file.path(local_benchmark_dir, "alevin-fry")
fs::dir_create(af_local_dir)
cellranger_local_dir <- file.path(local_benchmark_dir, "cellranger")
fs::dir_create(cellranger_local_dir)

# get list of directories containing results we want to copy over
# we only want the specific benchmarking runs that use splici, salign, and cr-like-em
af_dirs <- glue::glue("{alevin_s3_dir}/{benchmarking_run_ids}-Homo_sapiens.GRCh38.104.spliced_intron.txome-salign-cr-like-em")

# we only want alevin folder and any json files 
# copy to folder labeled with run id 
glue::glue(
  "aws s3 cp '{af_dirs}' '{af_local_dir}/{benchmarking_run_ids}' --exclude '*' --include 'alevin/*' --include '*.json' --recursive"
) |>
  purrr::walk(system)

# list of cellranger directories to copy over
cellranger_dirs <- glue::glue("{cellranger_s3_dir}/{benchmarking_run_ids}-GRCh38_104_cellranger_full-mRNA")

# we only need the filtered h5 file 
# copy each file to folder labeled with run id
glue::glue(
  "aws s3 cp '{cellranger_dirs}' '{cellranger_local_dir}/{benchmarking_run_ids}' --exclude '*' --include 'outs/filtered*.h5' --recursive"
) |> 
  purrr::walk(system)
