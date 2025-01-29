# This script calculates a TPM matrix for salmon quant files from all samples in a given project and exports an TSV file named `<project_id>-tpm.tsv`.
# Input quant files are expected to be stored as `project_id/sample_id/quant.sf` in the provided "input_dir" argument.
# Note that this script retains the original ensembl gene identifiers and does not perform gene symbol conversion

renv::load()
library(optparse)


# Parse options --------
option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    help = "Input directory with salmon quant files organized by project/sample"
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Output file to save TSV with TPM data."
  ),
  make_option(
    "--t2g_file",
    type = "character",
    default = here::here("s3_files", "reference_files", "Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv"),
    help = "Path to the t2g file used to read in quant.sf files."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs and paths -------
stopifnot(
  "An input directory must be provided to `input_dir`." = !is.null(opts$input_dir),
  "A path to an output file must be specified with `output_file`." = !is.null(opts$output_file),
  "The t2g file could not be found." = file.exists(opts$t2g_file)
)
stopifnot("Provided `input_dir` does not exist." = dir.exists(opts$input_dir))
fs::dir_create(dirname(opts$output_file))

# Read input helper files -----------
t2g_table <- readr::read_tsv(
  opts$t2g_file,
  col_names = c("transcript_id", "gene_id"),
  show_col_types = FALSE
)

# Find all quant.sf files -------------
quant_files <- list.files(
  path = opts$input_dir,
  pattern = "quant.sf",
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
tpm_df  <- txi_salmon$abundance |>
  as.data.frame() |>
  tibble::rownames_to_column("ensembl_id")


# Export to tsv -------
readr::write_tsv(tpm_df, opts$output_file)
