# This file contains functions used in prepare-scpca-portal-data.R to process consensus cell types.





#' Prepare and export TSV of consensus cell types
#'
#' @param sce SCE object to obtain cell types from 
#' @param output_tsv Path to output TSV file
#'
prepare_consensus_tsv <- function(sce, output_tsv){
  
  # Prepare and export the consensus cell types table from sce contents
  data.frame(
    project_id = metadata(sce)$project_id, 
    sample_id  = metadata(sce)$sample_id, 
    library_id = metadata(sce)$library_id, 
    sample_type = paste0(metadata(sce)$sample_type, collapse = ","),
    barcodes = colnames(sce), 
    # TODO TEMPORARY TO TEST THE CODE!!!!!!!!!!!
    consensus_annotation = sce$singler_celltype_annotation 
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
