# This script runs quanTIseq for samples in a given ScPCA project and exports 
#  a TSV of cell type proportions for each sample in the project.

renv::load()
library(optparse)
library(quantiseqr)

# Parse options --------
option_list <- list(
  make_option(
    "--project_id",
    type = "character",
    default = "SCPCP000001",
    help = "ScPCA project id to run quanTIseq on."
  ),
  make_option(
    "--input_dir",
    type = "character",
    default = here::here("analysis", "bulk-deconvolution", "data", "tpm"),
    help = "Input directory containing a project-specific TPM matrices saved as rds files named <project id>-tpm.rds."
  ),
  make_option(
    "--output_dir",
    type = "character",
    help = "Output directory to quanTIseq inference TSV."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs and define paths -------
stopifnot(
  "A project id in the format SCPCPXXXXXX must be provided to `project_id`." = !is.null(opts$project_id),
  "An output directory must be specified with `output_dir`." = !is.null(opts$output_dir)
)

tpm_file <- file.path(opts$input_dir, glue::glue("{opts$project_id}-tpm.rds"))
stopifnot("There is no TPM matrix for the given `project_id` in the provided `input_dir`." = file.exists(tpm_file))

fs::dir_create(opts$output_dir)
output_file <- file.path(opts$output_dir, glue::glue("{opts$project_id}-quantiseq.tsv"))


# Prepare input rds file -----------
tpm_matrix <- readr::read_rds(tpm_file)

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
  is_tumordata = FALSE, # as determined in ../exploratory-notebooks/quantiseq-tumor-genes.Rmd
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
readr::write_tsv(deconv_tidy_df, output_file)
