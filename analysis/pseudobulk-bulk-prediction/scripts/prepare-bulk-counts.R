# This script processes bulk counts expression into a list of data frames which can be used for analysis
# This will include some of the multiplexed samples, but these will be filtered out in next steps

renv::load()
library(optparse)
library(DESeq2)


# Parse options --------
option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    help = "Input directory containing all bulk TSV files, which will be found recursively"
  ),
  make_option(
    "--map_file",
    type = "character",
    help = "Path to TSV file mapping bulk sample and library ids"
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Path to output RDS containing normalized counts"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))


# Checks ------------------
stopifnot("Map file not found." = file.exists(opts$map_file))


# Find and read in bulk files --------
bulk_count_files <- list.files(
  path = here::here("analysis", "pseudobulk-bulk-prediction", "data", "scpca_data"), 
  full.names = TRUE,
  pattern = "_bulk_quant\\.tsv$", 
  recursive = TRUE
) 
stopifnot("No bulk count TSVs found" = length(bulk_count_files) > 0)
bulk_count_names <- stringr::str_split_i(basename(bulk_count_files), pattern = "_", i = 1)
names(bulk_count_files) <- bulk_count_names


# Read map to swap library to sample ids --------
sample_library_map <- readr::read_tsv(opts$map_file) |>
  dplyr::rename(
    sample_id = scpca_sample_id, 
    library_id = scpca_library_id
  )

# Make a list of data frames of bulk counts, normalized by DESeq2
bulk_counts_wide_df_list <- bulk_count_files |>
  purrr::map(
    \(counts_file) {
      
      df <- readr::read_tsv(counts_file, show_col_types = FALSE) 
      counts_mat <- as.matrix(df |> dplyr::select(-gene_id))
      rownames(counts_mat) <- df$gene_id
      counts_mat <- round(counts_mat) # deal with non-integers which DESeq2 will require
    
      bulk_counts_mat <- DESeqDataSetFromMatrix(
        countData = counts_mat,
        # these arguments don't matter for our purposes,
        # but DESeq2 requires them
        colData = data.frame(sample = colnames(counts_mat)),
        design = ~sample) |>
        estimateSizeFactors() |>
        rlog(blind = TRUE) |>
        assay() 
      
      bulk_counts_mat |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "gene_id") |>
        tidyr::pivot_longer(
          -gene_id, 
          names_to = "library_id", # library! 
          values_to = "bulk_counts"
        ) |>
        # swap to sample name
        dplyr::left_join(sample_library_map) |>
        dplyr::select(gene_id, sample_id, bulk_counts) |>
        tidyr::drop_na()
        
    }
  )

# Save to RDS ----------
readr::write_rds(bulk_counts_wide_df_list, opts$output_file)