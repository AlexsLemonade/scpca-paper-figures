# This file contains functions used in prepare-scpca-portal-data.R to process consensus cell types.





#' Prepare and export TSV of consensus cell types
#'
#' @param sce SCE object to obtain cell types from 
#' @param output_tsv Path to output TSV file
#'
prepare_consensus_tsv <- function(sce, output_tsv){
  
  # Define consensus cell types since they won't be present in SCEs which were not
  #  annotated with both SingleR and CellAssign
  if ("consensus_celltype_annotation" %in% names(colData(sce))) {
    # TODO: This column needs to be updated
    # See https://github.com/AlexsLemonade/scpca-paper-figures/issues/253
    consensus_celltypes <- sce$singler_celltype_annotation
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
  
  portal_metadata <- portal_metadata |>
    # Select columns of interest
    dplyr::select(
      scpca_project_id,
      scpca_sample_id,
      scpca_library_id,
      seq_unit, 
      technology
    ) |>
    # Group multiplexed sample ids back together
    dplyr::group_by(scpca_project_id, scpca_library_id, seq_unit, technology) |>
    dplyr::summarize(scpca_sample_id = paste(scpca_sample_id, collapse = ";"))
    
  
  # Define bulk and project metadata files
  bulk_files <- list.files(
    path = bulk_metadata_dir, 
    pattern = "*_bulk_metadata\\.tsv", 
    full.names = TRUE
  )
  project_files <- list.files(
    path = project_metadata_dir, 
    pattern = "*_single_cell_metadata\\.tsv", 
    full.names = TRUE
  ) |>
    # we need project id identifiers for this one
    purrr::set_names(\(x) {stringr::str_split_i(basename(x), "_", 1)})
  
  # Parse bulk metadata
  bulk_metadata <- bulk_files |>
    purrr::map(readr::read_tsv) |>
    purrr::list_rbind() |>
    dplyr::select(
      scpca_project_id = project_id, 
      scpca_sample_id = sample_id, 
      scpca_library_id = library_id, 
      seq_unit, 
      technology
    ) 
  
  # Parse project metadata to get rows for CITE-seq & cellhash libraries
  multimodal_metadata <- project_files |>
    purrr::imap(
      \(file, project_id) {
        metadata_df <- readr::read_tsv(file) 
        
        # Get information for cellhash and then CITE-seq projects
        if (project_id == "SCPCP000009") {
          metadata_df <- metadata_df |>
            dplyr::group_by(scpca_project_id, scpca_library_id, seq_unit, technology) |>
            dplyr::summarize(scpca_sample_id = paste(scpca_sample_id, collapse = ";")) |>
            dplyr::mutate(technology = "cellhash")
        } else {
          metadata_df <- metadata_df |>
            dplyr::filter(!is.na(adt_filtering_method)) |>
            dplyr::mutate(technology = "CITE-seq")
        } 
        
        metadata_df <- metadata_df |>
          dplyr::select(
            scpca_project_id,
            scpca_sample_id,
            scpca_library_id,
            seq_unit,
            technology
          )
        
        return(metadata_df)
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::ungroup()
  
  # Combine data frames and export
  portal_metadata |>
    dplyr::bind_rows(bulk_metadata) |>
    dplyr::bind_rows(multimodal_metadata) |>
    dplyr::arrange(scpca_project_id) |>
    readr::write_tsv(output_tsv)
}
