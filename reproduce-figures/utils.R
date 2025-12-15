# This file contains functions used in prepare-scpca-portal-data.R to process consensus cell types
#  and to create metadata files.

#' Prepare and export TSV of consensus cell types
#'
#' @param sce SCE object to obtain cell types from
#' @param output_tsv Path to output TSV file
#'
prepare_consensus_tsv <- function(sce, output_tsv){

  # Define consensus cell types since they won't be present in SCEs which were not
  #  annotated with both SingleR and CellAssign
  if ("consensus_celltype_annotation" %in% names(colData(sce))) {
    consensus_celltypes <- sce$consensus_celltype_annotation
  } else {
    consensus_celltypes <- NA
  }

  # Prepare and export the consensus cell types table from sce contents
  data.frame(
    project_id = metadata(sce)$project_id,
    sample_id  = metadata(sce)$sample_id,
    library_id = metadata(sce)$library_id,
    sample_type = paste0(metadata(sce)$sample_type, collapse = ","),
    barcodes = colnames(sce),
    consensus_annotation = consensus_celltypes
  ) |>
    # export
    readr::write_tsv(output_tsv)
}

#' Prepare and export TSV of marker gene expression
#'
#' This code is adapted from the consensus cell type module in OpenScPCA-nf:
#' https://github.com/AlexsLemonade/OpenScPCA-nf/blob/a97e88c411cabb1b41e1fce6a7735311fb435051/modules/cell-type-consensus/resources/usr/bin/assign-consensus-celltypes.R#L205
#'
#' @param sce SCE object to obtain gene expression data from
#' @param marker_genes Vector of marker genes to include expression for
#' @param output_tsv Path to output TSV file
#'
prepare_gene_expression_tsv <- function(sce, marker_genes, output_tsv) {

  expressed_markers <- rowData(sce) |>
    as.data.frame() |>
    dplyr::filter(detected > 0) |>
    dplyr::pull(gene_ids) |>
    intersect(marker_genes)

  # get logcounts from sce for expressed genes, unless the length is 0 in which case fill NAs
  if (length(expressed_markers) == 0) {
    gene_exp_df <- data.frame(
      barcodes = rownames(sce),
      library_id = metadata(sce)$library_id,
      ensembl_gene_id = NA,
      logcounts = NA
    )
  } else {
    gene_exp_df <- scuttle::makePerCellDF(
      sce,
      features = expressed_markers,
      assay.type = "logcounts",
      use.coldata = "barcodes",
      use.dimred = FALSE
    ) |>
      tidyr::pivot_longer(starts_with("ENSG"), names_to = "ensembl_gene_id", values_to = "logcounts") |>
      # add library id column
      dplyr::mutate(library_id = metadata(sce)$library_id, .before = 0)
  }

  # export the marker gene expression tsv
  readr::write_tsv(gene_exp_df, output_tsv)
}




#' Prepare and export sample metadata file
#'
#' @param portal_metadata Data frame of portal-wide metadata
#' @param output_tsv Path to output TSV file
#'
prepare_sample_metadata <- function(
    portal_metadata,
    output_tsv) {

  sample_metadata <- portal_metadata |>
    dplyr::select(
      scpca_project_id,
      scpca_sample_id,
      diagnosis,
      disease_timing,
      is_cell_line,
      is_xenograft) |>
    # remove duplicate rows, which occur when there are multiple libraries per sample
    dplyr::distinct()

  readr::write_tsv(sample_metadata, output_tsv)
}





#' Prepare and export library metadata file
#'
#' @param portal_metadata Data frame of portal-wide metadata
#' @param bulk_metadata_dir Directory with bulk metadata TSV files
#' @param project_metadata_dir Directory with project-specific metadata TSV files
#' @param output_tsv Path to output TSV file
#'
prepare_library_metadata <- function(
    portal_metadata,
    bulk_metadata_dir,
    project_metadata_dir,
    output_tsv) {

  # final columns to include in the exported library metadata tsv
  final_columns <- c(
    "scpca_project_id",
    "scpca_sample_id",
    "scpca_library_id",
    "seq_unit",
    "technology"
  )

  library_metadata <- portal_metadata |>
    # Group multiplexed sample ids back together; keep has_cellhash for later manipulation
    dplyr::group_by(scpca_project_id, scpca_library_id, seq_unit, technology, has_cellhash) |>
    dplyr::summarize(scpca_sample_id = paste(scpca_sample_id, collapse = ";")) |>
    dplyr::ungroup()

  # duplicate multiplexed rows so we have a row for cellhash technology
  cellhash_rows <- library_metadata |>
    dplyr::filter(has_cellhash) |>
    # in our code, we detect this technology with "cellhash" only, so the version isn't needed
    dplyr::mutate(technology = "cellhash")

  library_metadata <- library_metadata |>
    dplyr::bind_rows(cellhash_rows) |>
    dplyr::select(-has_cellhash)

  # Define bulk and project metadata files
  bulk_files <- list.files(
    path = bulk_metadata_dir,
    pattern = "*_bulk_metadata\\.tsv",
    full.names = TRUE
  )
  citeseq_files <- list.files(
    path = project_metadata_dir,
    pattern = "*_single_cell_metadata\\.tsv",
    full.names = TRUE
  )

  # Parse bulk metadata
  bulk_metadata <- bulk_files |>
    purrr::map(readr::read_tsv) |>
    purrr::list_rbind() |>
    dplyr::rename(
      scpca_project_id = project_id,
      scpca_sample_id = sample_id,
      scpca_library_id = library_id
    ) |>
    dplyr::select(all_of(final_columns))

  # Parse project metadata to get rows for CITE-seq
  cite_metadata <- citeseq_files |>
    purrr::map(
      \(file) {
        metadata_df <- readr::read_tsv(file) |>
          dplyr::filter(!is.na(adt_filtering_method)) |>
          # in our code, we detect this technology with "CITE" only, so the version isn't needed
          dplyr::mutate(technology = "CITEseq") |>
          dplyr::select(all_of(final_columns))

        return(metadata_df)
      }
    ) |>
    purrr::list_rbind()

  # Combine data frames and export
  library_metadata |>
    dplyr::bind_rows(bulk_metadata) |>
    dplyr::bind_rows(cite_metadata) |>
    readr::write_tsv(output_tsv)
}

