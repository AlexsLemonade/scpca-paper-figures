#!/usr/bin/env Rscript

# This script is used to download the consensus cell type results from OpenScPCA-nf
# All files will be saved in a folder named `s3_files/cell-type-consensus-results`

renv::load()

library(optparse)

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--profile"),
    type = "character",
    help = "name of AWS profile to use for copying files"
  ),
  make_option(
    opt_str = c("--openscpca_release"),
    type = "character",
    default = "2025-06-30",
    help = "OpenScPCA release date to download consensus cell type results from"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "Must provide an aws profile using --profile" = !is.null(opt$profile)
)

# read in the project whitelist
project_whitelist_file <- here::here("sample-info", "project-whitelist.txt")
project_ids <- readLines(project_whitelist_file)

# where openscpca results live
s3_results_bucket <- glue::glue("s3://openscpca-nf-workflow-results/{opt$openscpca_release}/cell-type-consensus"
all_s3_dirs <- glue::glue("{s3_results_bucket}/{project_ids}")

# specify local directories and create them
local_results_dir <- here::here("s3_files", "cell-type-consensus-results")
all_results_dir <- glue::glue("{local_results_dir}/{project_ids}")
fs::dir_create(all_results_dir)

# copy files with cell type results and gene expression tsvs for each sample
glue::glue(
  "aws s3 cp '{all_s3_dirs}' '{all_results_dir}' --exclude '*' --include '*_consensus-cell-types.tsv.gz' --include '*_marker-gene-expression.tsv.gz' --profile {opt$profile} --recursive"
) |>
  purrr::walk(system)
