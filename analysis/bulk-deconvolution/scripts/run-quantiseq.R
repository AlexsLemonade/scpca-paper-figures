# This script runs quanTIseq for samples in a given ScPCA project and exports
#  a TSV of cell type proportions for each sample in the project.
# Note that quanTIseq is deterministic and does not require a seed.

renv::load()
library(optparse)
library(quantiseqr)

# Parse options --------
option_list <- list(
  make_option(
    "--input_file",
    type = "character",
    help = "Input RDS file with a TPM matrix."
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Output TSV file to save quanTIseq inference to."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs and define paths -------
stopifnot(
  "A path to an input RDS file must be specified with `input_file`." = !is.null(opts$input_file),
  "A path to an output TSV must be specified with `output_file`." = !is.null(opts$output_file)
)

stopifnot("The provided `input_file` does not exist." = file.exists(opts$input_file))

fs::dir_create(dirname(opts$output_file))

# Prepare input rds file -----------
tpm_matrix <- readr::read_rds(opts$input_file)

# As determined in ../exploratory-notebooks/quantiseq-tumor-genes.Rmd, several gene names in the signature
#  are outdated. Here, we'll replace those gene symbols in our matrix row names with the gene symbols quanTIseq expects:
#  our data         quantiseq
#  PALM2AKAP2       AKAP2
#  TENT5C           FAM46C
#  GUCY1A1          GUCY1A3

new_rownames <- dplyr::case_match(
  rownames(tpm_matrix),
  "PALM2AKAP2" ~ "AKAP2",
  "TENT5C" ~ "FAM46C",
  "GUCY1A1" ~ "GUCY1A3",
  .default = rownames(tpm_matrix)
)
rownames(tpm_matrix) <- new_rownames

# Run quanTIseq -------------
deconv_df <- quantiseqr::run_quantiseq(
  expression_data = tpm_matrix,
  is_tumordata = TRUE, # as determined in ../exploratory-notebooks/quantiseq-tumor-genes.Rmd
  scale_mRNA = TRUE
)

# wrangle the output into a tidy format
deconv_tidy_df <- deconv_df |>
  # remove rownames
  tibble::as_tibble() |>
  tidyr::pivot_longer(
    -Sample,
    names_to = "quantiseq_celltype",
    values_to = "proportion"
  ) |>
  dplyr::rename(sample_id = Sample)

# Export to tsv -------
readr::write_tsv(deconv_tidy_df, opts$output_file)
