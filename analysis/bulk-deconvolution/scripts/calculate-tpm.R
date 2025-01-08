# This script calculates a TPM matrix for salmon quant files from all samples in a given project and exports an RDS file named `<project_id>-tpm.rds`.
# Input quant files are expected to be stored as `project_id/sample_id/library_id/quant.sf` in the provided "input_dir" argument.
# This file will have a column `gene_symbol` and a column for each sample.
# As such, this script also converts ensembl ids to gene symbols, which is needed as input to deconvolution methods.

renv::load()
library(optparse)


# Parse options --------
option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    help = "Input directory with salmon quant files organized by project/sample/library."
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Output file to save RDS with TPM matrix."
  ),
  make_option(
    "--t2g_file",
    type = "character",
    default = here::here("s3_files", "reference_files", "Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv"),
    help = "Path to the t2g file used to read in quant.sf files."
  ),
  make_option(
    "--id_map_file",
    type = "character",
    default = here::here("analysis", "bulk-deconvolution", "data", "reference", "ensembl_symbol.tsv"),
    help = "Path to TSV file with columns `ensembl_id` and `gene_symbol` which can be used to convert ids."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs and paths -------
stopifnot(
  "An input directory must be provided to `input_dir`." = !is.null(opts$input_dir),
  "A path to an output file must be specified with `output_file`." = !is.null(opts$output_file),
  "The t2g file could not be found." = file.exists(opts$t2g_file),
  "The id mapping file to use for id mapping could not be found." = file.exists(opts$id_map_file)
)
stopifnot("Provided `input_dir` does not exist." = dir.exists(opts$input_dir))
fs::dir_create(dirname(opts$output_file))

# Read input helper files -----------
t2g_table <- readr::read_tsv(
  opts$t2g_file,
  col_names = c("transcript_id", "gene_id"),
  show_col_types = FALSE
)
map_table <- readr::read_tsv(opts$id_map_file, show_col_types = FALSE)


# Find all quant.sf files -------------
quant_files <- list.files(
  path = opts$input_dir,
  recursive = TRUE
)
stopifnot("Could not find any quant.sf files for the specified project." = length(quant_files) > 0)

# Calculate TPM for each sample  ---------
sample_ids <- stringr::str_split_i(quant_files, pattern = "/", i = 1)
quant_paths <- setNames(file.path(opts$input_dir, quant_files), sample_ids)

txi_salmon <- tximport::tximport(
  quant_paths,
  type = "salmon",
  tx2gene = t2g_table
)
tpm_mat  <- txi_salmon$abundance

# Convert ensembl ids to gene symbols and sum TPM for duplicate symbols ------

# create a named vector for ensembl->symbol translation
symbol_map <- setNames(map_table$gene_symbol, map_table$ensembl_id)
# remove rows with no symbol
tpm_mat_pruned <- tpm_mat[!is.na(symbol_map[rownames(tpm_mat)]), ]
# sum across duplicated gene symbols
tpm_summed_mat <- rowsum(
  tpm_mat_pruned,
  group = symbol_map[rownames(tpm_mat_pruned)]
) # --> 39133 rows

# Export to rds -------
readr::write_rds(tpm_summed_mat, opts$output_file)
