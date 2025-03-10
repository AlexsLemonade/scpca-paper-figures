# This script prepares gene sets for input to over-representation analysis (ORA) for a given ScPCA project as follows:
# 
# 1. Identify all possible consensus cell types which _could_ have been called for the project
# 2. Identify the PanglaoDB cell types that contributed to those consensus labels
# 3. Using the scpca-nf PanglaoDB cell type reference for the given project, determine the genes in those PanglaoDB cell types
# 4. These (unique) genes define the gene set for a given consensus cell type to test in ORA
# 5. We include an indicator column `observed_in_singlecell` for whether the given gene set was observed in the project's consensus cell types 


renv::load()
library(optparse)


# Parse options --------
option_list <- list(
  make_option(
    "--project_id",
    type = "character",
    help = "Project being prepared"
  ),
  make_option(
    "--celltype_dir",
    type = "character",
    help = "Directory with consensus cell types for the given project"
  ),
  make_option(
    "--panglao_geneset_file",
    type = "character",
    help = "Path to gene set file used to cell type the given project"
  ), 
  make_option(
    "--output_file",
    type = "character",
    help = "Path to output TSV file to save gene sets."
  ), 
  make_option(
    "--consensus_reference_url",
    type = "character",
    default = "https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/refs/heads/main/analyses/cell-type-consensus/references/consensus-cell-type-reference.tsv",
    help = "URL with map between consensus cell types and PanglaoDB gene set labels"
  ), 
  make_option(
    "--library_metadata_file",
    type = "character",
    default = here::here("s3_files", "scpca-library-metadata.tsv"),
    help = "Path to ScPCA library metadata file."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Checks---------
stopifnot(
  "Panglao gene set file not found" = file.exists(opts$panglao_geneset_file), 
  "Library metadata file not found" = file.exists(opts$library_metadata_file)
)

# Prepare input files  -----------------
consensus_ref_df <- readr::read_tsv(opts$consensus_reference_url)

panglao_df <- readr::read_tsv(opts$panglao_geneset_file) |>
  tidyr::pivot_longer(
    -ensembl_id, 
    names_to = "panglao_celltype", 
    values_to = "gene_present"
  )

celltype_files <- list.files(
  opts$celltype_dir, 
  pattern = "_processed_consensus-cell-types\\.tsv\\.gz$", 
  full.names = TRUE, 
  recursive = TRUE
) |>
  purrr::set_names(
    \(x) {stringr::str_split_i(basename(x), pattern = "_", i = 1)}
  )

library_metadata <- readr::read_tsv(opts$library_metadata_file)

# determine libraries to keep from the consensus cell types
# we need to take this step because some projects don't use all their single-cell libraries
bulk_samples <- library_metadata |>
  dplyr::filter(
    scpca_project_id == opts$project_id,
    seq_unit == "bulk"
  ) |>
  dplyr::pull(scpca_sample_id)

keep_libraries <- library_metadata |>
  dplyr::filter(
    scpca_sample_id %in% bulk_samples,
    !stringr::str_detect(scpca_sample_id, ";"),
    seq_unit %in% c("cell", "nucleus")
  ) |>
  dplyr::pull(scpca_library_id) |>
  unique()

# Find all the _observed_ consensus cell type files for this project
# we can use these to add an indicator to the final gene set list
observed_celltypes <- celltype_files |>
  # only keep relevant libraries
  purrr::keep_at(
    \(library_id) library_id %in% keep_libraries
  ) |>
  purrr::map(
    \(f) {
      readr::read_tsv(f) |>
        dplyr::pull(consensus_annotation) 
    }) |>
  purrr::reduce(c) |>
  unique()
observed_celltypes <- observed_celltypes[observed_celltypes != "Unknown"]


# Determine consensus gene sets -------------------

# Find all panglao cell types in the Panglao reference
panglao_celltypes <- unique(panglao_df$panglao_celltype)
panglao_celltypes <- panglao_celltypes[panglao_celltypes != "other"]


# Determine the consensus cell types they contribute to
consensus_cell_types <- consensus_ref_df |>
  dplyr::filter(original_panglao_name %in% panglao_celltypes) |>
  dplyr::select(consensus_annotation, original_panglao_name) |>
  dplyr::distinct()

# Create the consensus gene sets 
consensus_genesets_df <- split(consensus_cell_types, consensus_cell_types$consensus_annotation) |>
  purrr::map(
    \(df) {
      panglao_df |>
        dplyr::filter(
          panglao_celltype %in% df$original_panglao_name, 
          gene_present == 1
        ) |>
        dplyr::select(ensembl_id) |>
        dplyr::distinct()
    }
  ) |>
  # Convert to TSV for human-readable export
  purrr::list_rbind(names_to = "cell_type_name") |>
  # Now, add an indicator column for whether the gene set (cell type) was observed
  dplyr::mutate(observed_in_singlecell = cell_type_name %in% observed_celltypes)



# Export tsv file ------------------------------
readr::write_tsv(consensus_genesets_df, opts$output_file)
