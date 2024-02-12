# This script is used to download any data needed to generate figures 
# All files will be saved in a folder named `s3_files`

renv::load()

# sync results files for SCPCS000001 (for QC plots + cell type plots) ----------

qc_s3_files <- "s3://nextflow-ccdl-results/scpca-prod/results/SCPCP000001/SCPCS000001"
qc_local_files <- here::here("s3_files", "SCPCS000001")
fs::dir_create(qc_local_files)

sync_call <- glue::glue("aws s3 sync '{qc_s3_files}' '{qc_local_files}' --exclude '*' --include '*.rds' --exact-timestamps")
system(sync_call)

# sync results for SCPCS000216 (for ADT plots) ---------------------------------

adt_s3_files <- "s3://nextflow-ccdl-results/scpca-prod/results/SCPCP000007/SCPCS000216"
adt_local_files <- here::here("s3_files", "SCPCS000216")
fs::dir_create(adt_local_files)

sync_call <- glue::glue("aws s3 sync '{adt_s3_files}' '{adt_local_files}' --exclude '*' --include '*.rds' --exact-timestamps")
system(sync_call)

# sync results for benchmarking libraries --------------------------------------

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
  "aws s3 sync '{af_dirs}' '{af_local_dir}/{benchmarking_run_ids}' --exclude '*' --include 'alevin/*' --include '*.json' --exact-timestamps"
) |>
  purrr::walk(system)

# list of cellranger directories to copy over
cellranger_dirs <- glue::glue("{cellranger_s3_dir}/{benchmarking_run_ids}-GRCh38_104_cellranger_full-mRNA")

# we only need the filtered h5 file 
# copy each file to folder labeled with run id
glue::glue(
  "aws s3 sync '{cellranger_dirs}' '{cellranger_local_dir}/{benchmarking_run_ids}' --exclude '*' --include 'outs/filtered*.h5' --exact-timestamps"
) |> 
  purrr::walk(system)

# sync results files for SCPCL000498 (for submitter cell type heatmap) ---------
celltype_s3_files <- "s3://nextflow-ccdl-results/scpca-prod/results/SCPCP000005/SCPCS000251"
celltype_local_files <- here::here("s3_files", "SCPCS000251")
fs::dir_create(celltype_local_files)

sync_call <- glue::glue("aws s3 sync '{celltype_s3_files}' '{celltype_local_files}' --exclude '*' --include 'SCPCL000498_processed.rds' --exact-timestamps")
system(sync_call)

# sync merged object results for SCPCP000003 -----------------------------------
merged_s3_dir <- "s3://nextflow-ccdl-results/scpca/processed/results/merged/SCPCP000003"
merged_local_dir <- here::here("s3_files", "SCPCP000003")
fs::dir_create(merged_local_dir)

merged_file_name <- "SCPCP000003_merged.rds"

sync_call <- glue::glue("aws s3 sync '{merged_s3_dir}' '{merged_local_dir}' --exclude '*' --include '{merged_file_name}' --exact-timestamps")
system(sync_call)
