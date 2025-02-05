# This script holds helper functions used by notebooks in the `model-notebooks` directory


#' Prepare data for analysis
#'
#' This function reads and combines bulk and pseudobulk data into a single
#' data frame for a given project
#'
#' @param bulk_counts_file Path to RDS file with bulk counts for all projects
#' @param pseudobulk_file Path to pseudobulk TSV file for the project of interest
#' @param project_id Project id of interest, used to filter the bulk_counts
#' @param exclude_samples Optional vector of samples to remove from the final data frame
#'
#' @returns
#' @export
#'
#' @examples
prepare_data <- function(
    bulk_counts_file,
    pseudobulk_file,
    project_id,
    exclude_samples = c()
) {


  bulk_counts_df <- readr::read_rds(bulk_counts_file) |>
    purrr::pluck(project_id) |>
    # make names consistent with pseudobulk df
    dplyr::rename(
      bulk = bulk_counts,
      ensembl_id = gene_id
    )

  pseudo_df <- readr::read_tsv(pseudobulk_file, show_col_types = FALSE) |>
    dplyr::filter(expression_type == "pseudobulk_deseq") |>
    dplyr::select(ensembl_id, sample_id, pseudobulk = expression)

  # combine bulk and pseudobulk into a single data frame
  expr_df <- dplyr::full_join(
    pseudo_df,
    bulk_counts_df,
    by = c("ensembl_id", "sample_id")
  ) |>
    dplyr::filter(!(sample_id %in% exclude_samples))

  return(expr_df)
}
