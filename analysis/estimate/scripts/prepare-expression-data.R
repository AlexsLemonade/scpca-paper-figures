# This script prepares data for a given project for analysis, including:
# - Bulk counts
# - Bulk TPM
# - Pseudobulk counts
# For each measure, we create a list of data frames for each project in the wide format expected by `tidyestimate`
# No normalization or transformation is performed since ESTIMATE bases analysis on ranks


renv::load()
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})


# Parse options --------
option_list <- list(
  make_option(
    "--project_id",
    type = "character",
    help = "ScPCA project ID to prepare data for"
  ),
  make_option(
    "--scpca_dir",
    type = "character",
    help = "Input project directory with ScPCA data files to find recursively"
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Path to RDS file with list of data frames of expression quantites"
  ),
  make_option(
    "--library_metadata_file",
    type = "character",
    default = here::here("s3_files", "scpca-library-metadata.tsv"),
    help = "Path to ScPCA library metadata file"
  ),
  make_option(
    "--ensembl_symbol_map_file",
    type = "character",
    default = here::here("analysis", "estimate", "data", "ensembl-symbol-map.tsv"),
    help = "Path to TSV file with map between ensembl ids and gene symbols"
  ),
  make_option(
    "--t2g_file",
    type = "character",
    default = here::here("s3_files", "reference_files", "Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv"),
    help = "Path to the t2g file used to read in quant.sf files."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Functions --------------------------

# Helper function to convert ensembl ids to gene symbols
# It also sums expression, per sample, of duplicated gene symbols
convert_ids <- function(df, id_map_df) {
  
  summed_df <- df |>
    #### use gene symbols instead of ensembl
    dplyr::inner_join(id_map_df, by = "ensembl_id") |>
    tidyr::drop_na(gene_symbol) |>
    #### sum duplicates
    tidyr::pivot_longer(
      contains("SCPCS"), 
      names_to = "sample_id", 
      values_to = "expr") |>
    dplyr::group_by(gene_symbol, sample_id) |>
    dplyr::summarize(total = sum(expr)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      names_from = sample_id, 
      values_from = total
    ) |>
    #### order columns as expected by tidyestimate
    dplyr::select(gene_symbol, contains("SCPCS"))
  
  stopifnot("Duplicate gene ids not successfully handled" = !any(duplicated(summed_df$gene_symbol)))
  
  return(summed_df)
}

# Checks and file setup ------------------
stopifnot(
  "The ScPCA data directory must be provided to `scpca_dir`." = !is.null(opts$scpca_dir),
  "A path to an output RDS file must be specified with `output_file`." = !is.null(opts$output_file), 
  "Library metadata file not found." = file.exists(opts$library_metadata_file),
  "ID map file not found." = file.exists(opts$ensembl_symbol_map_file),
  "t2g file not found." = file.exists(opts$t2g_file)
)

fs::dir_create(dirname(opts$output_file))

# Read ensembl <-> symbol map file
id_map_df <- readr::read_tsv(opts$ensembl_symbol_map_file, show_col_types = FALSE)



# Prepare TPM counts data frame ----------------

t2g_table <- readr::read_tsv(
  opts$t2g_file,
  col_names = c("transcript_id", "gene_id"),
  show_col_types = FALSE
)

# Find all quant.sf files 
quant_files <- list.files(
  path = opts$scpca_dir,
  pattern = "quant.sf",
  recursive = TRUE
)
stopifnot("Could not find any quant.sf files for the specified project." = length(quant_files) > 0)

# Calculate TPM for each sample 
sample_ids <- stringr::str_split_i(quant_files, pattern = "/", i = 1)
quant_paths <- setNames(file.path(opts$scpca_dir, quant_files), sample_ids)

tpm_df <- tximport::tximport(
  quant_paths,
  type = "salmon",
  tx2gene = t2g_table
) |>
  purrr::pluck("abundance") |>
  as.data.frame() |>
  tibble::rownames_to_column("ensembl_id") |>
  convert_ids(id_map_df)

tpm_samples <- colnames(tpm_df[,-1]) # used to prepare bulk counts

# Prepare bulk counts data frame ----------------

bulk_file <- file.path(
  opts$scpca_dir, 
  glue::glue("{opts$project_id}_bulk_quant.tsv")
)
raw_bulk_df <- readr::read_tsv(bulk_file, show_col_types = FALSE) 

# read library metadata to swap library to sample ids
sample_library_df <- readr::read_tsv(opts$library_metadata_file, show_col_types = FALSE) |>
  dplyr::filter(
    scpca_project_id == opts$project_id,
    seq_unit == "bulk"
  ) |> 
  dplyr::select(scpca_sample_id, scpca_library_id)

bulk_df <- raw_bulk_df |>
  dplyr::rename(ensembl_id = gene_id) |>
  # replace library with sample ids
  tidyr::pivot_longer(contains("SCPCL")) |>
  dplyr::left_join(sample_library_df, by = c("name" = "scpca_library_id")) |>
  dplyr::select(-name) |>
  # subset to only bulk samples in the TPM df, since bulk raw counts may have more columns than we want
  dplyr::filter(scpca_sample_id %in% tpm_samples) |>
  tidyr::pivot_wider(
    names_from = scpca_sample_id, 
    values_from = value
  ) |>
  convert_ids(id_map_df)

# Prepare pseudobulk counts data frame ----------------

# read all sces
rds_files <- list.files(
  path = opts$scpca_dir,
  pattern = "_processed\\.rds$",
  recursive = TRUE
)
stopifnot(
  "Could not find any _processed.rds files in the provided `input_dir`." = length(rds_files) > 0
)
sample_ids <- stringr::str_split_i(rds_files, pattern = "/", i = 1)
rds_files <- file.path(opts$scpca_dir, rds_files) # update to contain full path
names(rds_files) <- sample_ids
sce_list <- purrr::map(rds_files, readr::read_rds)

# extract and sum the raw counts
pseudo_df <- sce_list |>
  purrr::map(counts) |>
  purrr::map(rowSums) |>
  do.call(cbind, args = _ ) |>
  as.data.frame() |>
  tibble::rownames_to_column("ensembl_id") |>
  convert_ids(id_map_df)



# Export RDS -----------------

expr_list <- list(
  "bulk_counts" = bulk_df, 
  "bulk_tpm" = tpm_df, 
  "pseudobulk" = pseudo_df
)
# check that all data frames are compatible:

# _bulk_ genes should be equal
stopifnot("Bulk data frames have different genes" = nrow(bulk_df) == nrow(tpm_df))

# samples should match all around
stopifnot(
  "Data frames have different samples" = purrr::every(expr_list, \(x) all(colnames(x) == colnames(bulk_df)))
)

readr::write_rds(expr_list, opts$output_file)
