# This script runs EPIC, using both of its built-in gene signatures sets,
# for all bulk samples in a given ScPCA project.
# It exports a TSV of cell type proportions for each sample in the project
# for both the TRef and BRef reference in EPIC.
# Note that EPIC does not require a seed.

renv::load()
library(optparse)
library(EPIC)

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
    help = "Output TSV file to save EPIC inferences to."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Define functions --------------------

# This is a helper function to format the cell fraction matrix that EPIC infers
# into a long-format data frame with additional columns for:
# - the reference (`ref_name`)
# - the `mRNAProportions` fractions
# - whether the model for the given sample converged
format_epic_output <- function(epic_output, ref_name) {

  # First, make a data frame of sample convergence
  convergence_df <- epic_output |>
    purrr::pluck("fit.gof") |>
    tibble::rownames_to_column(var = "sample_id") |>
    # the convergeCode is like bash: 1 is fail, 0 is success
    dplyr::mutate(converged = !as.logical(convergeCode)) |>
    dplyr::select(sample_id, converged)
  
  # Extract `cellFractions` into long data frame
  epic_df <- epic_output |>
    purrr::pluck("cellFractions") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "sample_id") |>
    tidyr::pivot_longer(
      -sample_id,
      names_to = "epic_celltype",
      values_to = "fraction"
    )
  
  # Add column with mRNAProportions
  epic_output |>
    purrr::pluck("mRNAProportions") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "sample_id") |>
    # name columns to facilitate joining
    tidyr::pivot_longer(
      -sample_id,
      names_to = "epic_celltype", # this column will be used to join
      values_to = "mRNAProportions" # this column will actually contain the values
    ) |>
    # join back into epic_df
    dplyr::right_join(epic_df, by = c("sample_id", "epic_celltype")) |>
    # arrange columns
    dplyr::select(sample_id, epic_celltype, fraction, mRNAProportions)

}

# Check inputs and define paths -------
stopifnot(
  "A path to an input RDS file must be specified with `input_file`." = !is.null(opts$input_file),
  "A path to an output TSV must be specified with `output_tsv_file`." = !is.null(opts$output_file)
)

stopifnot("The provided `input_file` does not exist." = file.exists(opts$input_file))

fs::dir_create(dirname(opts$output_file))

# Prepare input rds file -----------
tpm_matrix <- readr::read_rds(opts$input_file)

# As determined in ../exploratory-notebooks/epic-signature-genes.Rmd, several
#  gene names in the signatures are outdated. Here, we'll replace those gene symbols
#  in our matrix row names with the gene symbols EPIC expects:
#  our data     EPIC
# --------------------------
#  DIPK2B       CXorf36
#  ADA2         CECR1
#  H2BC5        HIST1H2BC
#  H3C4         HIST1H3D
new_rownames <- dplyr::case_match(
  rownames(tpm_matrix),
  "DIPK2B" ~ "CXorf36",
  "ADA2" ~ "CECR1",
  "H2BC5" ~ "HIST1H2BC",
  "H3C4" ~ "HIST1H3D",
  .default = rownames(tpm_matrix)
)
rownames(tpm_matrix) <- new_rownames

# Run EPIC with both references  -------------

ref_list <- c("TRef", "BRef")
names(ref_list) <- ref_list
epic_list <- ref_list |>
  purrr::map(\(ref){
    EPIC(bulk = tpm_matrix, reference = ref)
  })

# Format and combine fractions into single data frame
epic_df <- epic_list |>
  purrr::imap(format_epic_output) |>
  dplyr::bind_rows()

# Export results-------
readr::write_tsv(epic_df, opts$output_file)

