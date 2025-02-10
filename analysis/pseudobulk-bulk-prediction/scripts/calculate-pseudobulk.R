# This script calculates two flavors of pseudobulk normalized counts for a given project.
# The script exports an RDS file with a list of two matrices:
#  1. `pseudobulk_deseq`: Pseudobulk calculated by summing raw counts and normalizing them with DESeq2
#  2. `pseudobulk_log_counts`: Pseudobulk calculated by summing counts and taking their log2
#
# The script also exports a TSV of the percent of samples each gene is expressed in

renv::load()
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(SingleCellExperiment)
})


# Parse options --------
option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    help = "Input directory containing all _processed.rds files, which will be found recursively"
  ),
  make_option(
    "--output_pseudobulk_file",
    type = "character",
    help = "Path to output TSV to save pseudobulk calculations"
  ),
  make_option(
    "--output_percent_expressed_file",
    type = "character",
    help = "Path to output TSV to save percent of samples each gene is expressed in, based on raw counts"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs and paths -------
stopifnot(
  "An input directory must be provided to `input_dir`." = !is.null(opts$input_dir),
  "A path to an output TSV file to save pseudobulk counts must be specified with `output_pseudoulk_file`." = !is.null(opts$output_pseudobulk_file),
  "A path to an output TSV file to save percent of single-cell samples genes are expressed in must be specified with `output_percent_expressed_file`." = !is.null(opts$output_percent_expressed_file)
)
fs::dir_create(dirname(opts$output_pseudobulk_file))
fs::dir_create(dirname(opts$output_percent_expressed_file))


# Read in all SCEs----------------
rds_files <- list.files(
  path = opts$input_dir,
  pattern = "_processed\\.rds$",
  recursive = TRUE
)
stopifnot(
  "Could not find any _processed.rds files in the provided `input_dir`." = length(rds_files) > 0
)

sample_ids <- stringr::str_split_i(rds_files, pattern = "/", i = 1)
rds_files <- file.path(opts$input_dir, rds_files) # update to contain full path
names(rds_files) <- sample_ids

sce_list <- purrr::map(rds_files, readr::read_rds)

# extract and sum the raw counts
pseudo_raw_counts <- sce_list |>
  purrr::map(counts) |>
  purrr::map(rowSums) |>
  do.call(cbind, args = _ )


# Export table of percent expressed based on raw counts alone
rowMeans(pseudo_raw_counts > 0) |>
  tibble::as_tibble(rownames = "ensembl_id") |>
  dplyr::rename(percent_samples_expressed = value) |>
  readr::write_tsv(opts$output_percent_expressed_file)


# Approach 1: Normalize with DESeq 2 ----------------
pseudo_deseq <- DESeqDataSetFromMatrix(
  countData = pseudo_raw_counts,
  # these arguments don't matter for our purposes,
  # but DESeq2 requires them
  colData = data.frame(sample = sample_ids),
  design = ~sample) |>
  estimateSizeFactors() |>
  rlog(blind = TRUE) |>
  assay() 


# Approach 2: Sum counts and log1p (but base 2) directly ----------------
pseudo_log_counts <- log1p(pseudo_raw_counts)/log(2)

# Combine into a single long data frame ---------------------

# Define a helper function for converting each matrix to a long data frame
make_long <- function(pseudo_mat, label) {
  pseudo_mat |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "ensembl_id") |>
    tidyr::pivot_longer(
      -ensembl_id, 
      names_to = "sample_id", 
      values_to = "expression"
    ) |>
    dplyr::mutate(expression_type = label)
}

pseudo_df <- dplyr::bind_rows(
  make_long(pseudo_deseq, "pseudobulk_deseq"), 
  make_long(pseudo_log_counts, "pseudobulk_log_counts")
)

# Export ------------------
readr::write_tsv(pseudo_df, opts$output_pseudobulk_file)
