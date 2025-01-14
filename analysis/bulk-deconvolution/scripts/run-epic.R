# This script runs EPIC, using both of its built-in gene signatures sets,
# for all bulk samples in a given ScPCA project.
# It a exports a TSV of cell type proportions for each sample in the project 
#  and an RDS file with a list of the full EPIC objects, named by reference

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
    "--output_tsv_file",
    type = "character",
    help = "Output TSV file to save EPIC inferences to."
  ), 
  make_option(
    "--output_epic_object_file",
    type = "character",
    help = "Output RDS file to save the full EPIC objects to."
  )   
)
opts <- parse_args(OptionParser(option_list = option_list))


# Define functions --------------------

# This is a helper function to format the cell fraction matrix that EPIC infers
# into a long-format data frame with columns indicating the reference and whether
# the sample converged
format_epic_fractions <- function(epic_output, ref_name) {
  
  # First, make a data frame of sample convergence
  convergence_df <- epic_output |>
    purrr::pluck("fit.gof") |>
    tibble::rownames_to_column(var = "sample_id") |>
    # the convergeCode is like bash: 1 is fail, 0 is success
    dplyr::mutate(converged = !as.logical(convergeCode)) |>
    dplyr::select(sample_id, converged)
  
  epic_output |>
    purrr::pluck("cellFractions") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "sample_id") |>
    tidyr::pivot_longer(
      -sample_id,
      names_to = "cell_type", 
      values_to = "fraction"
    ) |>
    dplyr::mutate(reference = ref_name) |>
    dplyr::full_join(convergence_df, by = "sample_id")
}


# Check inputs and define paths -------
stopifnot(
  "A path to an input RDS file must be specified with `input_file`." = !is.null(opts$input_file),
  "A path to an output TSV must be specified with `output_tsv_file`." = !is.null(opts$output_tsv_file),
  "A path to an output RDS file must be specified with `output_epic_object_file`." = !is.null(opts$output_epic_object_file)
)

stopifnot("The provided `input_file` does not exist." = file.exists(opts$input_file))

fs::dir_create(dirname(opts$output_tsv_file))
fs::dir_create(dirname(opts$output_epic_object_file))

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

epic_tref <- EPIC(
  bulk = tpm_matrix, 
  reference = "TRef"
) 

epic_bref <- EPIC(
  bulk = tpm_matrix, 
  reference = "BRef"
) 

epic_list <- list(
  "Tref" = epic_tref, 
  "Bref" = epic_bref
)


# Format and combine fractions into single data frame
epic_df <- epic_list |>
  purrr::imap(format_epic_fractions) |>
  dplyr::bind_rows()

# Export results-------

# fractions
readr::write_tsv(epic_df, opts$output_tsv_file)

# full result
readr::write_rds(epic_list, opts$output_epic_object_file)