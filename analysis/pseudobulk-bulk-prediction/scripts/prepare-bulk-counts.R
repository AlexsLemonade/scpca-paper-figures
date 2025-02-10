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
    "--output_counts_file",
    type = "character",
    help = "Path to output RDS containing normalized counts"
  ),
  make_option(
    "--output_percent_expressed_file",
    type = "character",
    help = "Path to output TSV to save, per project, the percent of samples each gene is expressed in, based on raw counts"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))


# Checks ------------------
stopifnot(
  "An input directory must be provided to `input_dir`." = !is.null(opts$input_dir),
  "Map file not found." = file.exists(opts$map_file),
  "A path to an output TSV file to save counts must be specified with `output_counts_file`." = !is.null(opts$output_counts_file),
  "A path to an output TSV file to save percent per project of single-cell samples genes are expressed in must be specified with `output_percent_expressed_file`." = !is.null(opts$output_percent_expressed_file)
)
fs::dir_create(dirname(opts$output_counts_file))
fs::dir_create(dirname(opts$output_percent_expressed_file))


# Find and read in bulk files --------
bulk_count_files <- list.files(
  path = opts$input_dir, 
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
 bulk_dat <- bulk_count_files |>
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
      
      # data frame of bulk counts 
      bulk_counts_df <- bulk_counts_mat |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "ensembl_id") |>
        tidyr::pivot_longer(
          -ensembl_id, 
          names_to = "library_id", # library! 
          values_to = "bulk_counts"
        ) |>
        # swap to sample name
        dplyr::left_join(sample_library_map) |>
        dplyr::select(ensembl_id, sample_id, bulk_counts) |>
        tidyr::drop_na()
      
      # data frame of percent of samples genes are expressed in
      percent_expr_df <- rowMeans(counts_mat > 0) |>
        tibble::as_tibble(rownames = "ensembl_id") |>
        dplyr::rename(percent_samples_expressed = value) 
      
      return(
        list(bulk_counts_df,  percent_expr_df)
      )
    }
  ) |> 
   purrr::transpose()

# Export -----------
 
# Save counts to RDS
readr::write_rds(bulk_dat[[1]], opts$output_counts_file)

# Save percent expressed to TSV
bulk_dat[[2]] |> 
  purrr::list_rbind(names_to = "project_id") |>
  readr::write_rds(opts$output_percent_expressed_file)
