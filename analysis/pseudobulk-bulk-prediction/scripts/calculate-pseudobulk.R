# This script calculates two flavors of pseudobulk normalized counts for a given project.
# The script exports an RDS file with a list of two matrices:
#  1. `pseudobulk_deseq`: Pseudobulk calculated by summing raw counts and normalizing them with DESeq2
#  2. `pseudobulk_logcounts`: Pseudobulk calculated by summing logcounts directly
#
# Currently, we also discard genes for which

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
    "--output_file",
    type = "character",
    help = "Path to output RDS to save pseudobulk matrices"
  ),
  make_option(
    "--min_count_per_gene",
    type = "integer",
    default = 10, # TODO: Do we just to get rid of rows that are all 0s, and remove this option entirely?
    help = "Threshold raw expression level to retain genes. If the pseudobulk count is less than this value, we discard the gene."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs and paths -------
stopifnot(
  "An input directory must be provided to `input_dir`." = !is.null(opts$input_dir),
  "A path to an output file must be specified with `output_file`." = !is.null(opts$output_file)
)
fs::dir_create(dirname(opts$output_file))


# Read in all SCEs----------------
rds_files <- list.files(
  path = opts$input_dir,
  pattern = "*_processed.rds",
  recursive = TRUE
)
stopifnot(
  "Could not find any _processed.rds files in the provided `input_dir`." = length(rds_files) > 0
)

sample_ids <- stringr::str_split_i(rds_files, pattern = "/", i = 1)
rds_files <- file.path(opts$input_dir, rds_files) # update to contain full path
names(rds_files) <- sample_ids

sce_list <- purrr::map(rds_files, readr::read_rds)

pseudo_raw_counts <- sce_list |>
  counts() |>
  purrr::map(DelayedArray::rowSums) |>
  # reduce seems to make the names go away, which is sad
  purrr::reduce(cbind)

colnames(pseudo_raw_counts) <- sample_ids

# We'll also identify rows to keep based on opts$min_count_per_gene
keep_rows <- DelayedArray::rowSums(pseudo_raw_counts >= opts$min_count_per_gene) > 0
pseudo_raw_counts <- pseudo_raw_counts[keep_rows, ]


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


# Approach 2: Sum logcounts directly ----------------
pseudo_logcounts <- sce_list |>
  purrr::map(
    \(sce) {
      # subset to the `keep_rows` we previously identified from the raw counts
      DelayedArray::rowSums(logcounts(sce)[keep_rows, ])
    }
  ) |>
  purrr::reduce(cbind)
colnames(pseudo_logcounts) <- sample_ids


# Export ------------------

# Confirm we have the same genes
stopifnot(
  setequal(rownames(pseudo_deseq), rownames(pseudo_logcounts))
)

pseudo_list <- list(
  pseudobulk_deseq = pseudo_deseq,
  pseudobulk_logcounts = pseudo_logcounts
)
readr::write_rds(pseudo_list, opts$output_file)