# Helper functions for creating cell type plots 
# These functions are taken directly from `scpca-nf/templates/celltypes_qc.rmd`,
# `scpca-nf/templates/celltypes_supplemental_report.rmd`, and `scpca-nf/templates/utils/celltype_functions.Rmd`

library(SingleCellExperiment)

# Set up cell type dataframe ---------------------------------------------------

#' Create `celltype_df` data frame for use in cell type QC reports
#'
#' @param processed_sce The processed sce object with cell type annotations in colData
#'
#' @return `celltype_df` with column of cell types, as factors, for each annotation method
create_celltype_df <- function(processed_sce) {
  # only incorporate UMAP coordinates if present
  if ("UMAP" %in% reducedDimNames(processed_sce)) {
    celltype_df <- processed_sce |>
      scuttle::makePerCellDF(use.dimred = "UMAP") |>
      # rename UMAP columns as needed to remove potential period added by `scuttle::makePerCellDF`
      dplyr::rename_with(
        \(x) stringr::str_replace(x, "^UMAP\\.", "UMAP"),
        starts_with("UMAP")
      )
    # otherwise just grab the colData
  } else {
    celltype_df <- colData(processed_sce) |>
      as.data.frame()
  }
  
  celltype_df <- celltype_df |>
    # only keep columns of interest
    dplyr::select(
      barcodes,
      # account for potentially missing columns
      contains("cluster"),
      contains("UMAP"),
      contains("singler"),
      contains("cellassign"),
      contains("submitter")
    )
  
  if ("submitter_celltype_annotation" %in% names(celltype_df)) {
    celltype_df <- prepare_submitter_annotation_values(celltype_df)
  }
  
  if ("singler_celltype_annotation" %in% names(celltype_df)) {
    celltype_df <- prepare_automated_annotation_values(
      celltype_df,
      singler_celltype_annotation
    )
  }
  if ("cellassign_celltype_annotation" %in% names(celltype_df)) {
    celltype_df <- prepare_automated_annotation_values(
      celltype_df,
      cellassign_celltype_annotation
    )
  }
  
  return(celltype_df)
}


#' Prepare and reformat cell type automated annotation values for use in QC reports
#'  Unknown cell types are updated with the label "Unknown cell type", and
#'  cell types are ordered in order of descending frequency, but with
#'  "Unknown cell type" as the last level
#'
#' @param df The data frame containing cell type annotations, one row per cell
#' @param annotation_column The column (written plainly, not a string) containing annotations to reformat
#'
#' @return Updated data frame with the `annotation_column` reformatted
prepare_automated_annotation_values <- function(
    df,
    annotation_column) {
  unknown_string <- "Unknown cell type"
  
  df |>
    dplyr::mutate(
      {{ annotation_column }} := dplyr::case_when(
        # singler condition
        is.na({{ annotation_column }}) ~ unknown_string,
        # cellassign conditon
        {{ annotation_column }} == "other" ~ unknown_string,
        # otherwise, keep it
        .default = {{ annotation_column }}
      ) |>
        # order column by number of cells
        forcats::fct_infreq() |>
        # make "Unknown cell type" the last level
        forcats::fct_relevel(unknown_string, after = Inf)
    )
}

# Jaccard Calculation ----------------------------------------------------------

#' Function to calculate Jaccard similarity on two vectors
#'
#' @param vec1 First vector
#' @param vec2 Second vector
#'
#' @return Jaccard similarity between the vectors
jaccard <- function(vec1, vec2) {
  length(intersect(vec1, vec2)) / length(union(vec1, vec2))
}


# Wrapper function to calculate jaccard similarity matrix for two categorical variables
#'
#' @param celltype_df The celltype_df data frame which must contain these columns:
#'   `colname1`, `colname2`, and `barcodes`
#' @param colname1 Column name, as a string, of first categorical variable of interest
#' @param colname2 Column name, as a string, of second categorical variable of interest
#'
#' @return Jaccard similarity matrix for the two columns. `colname1` values will
#'   be row names and `colname2` values will be column names in the final matrix
make_jaccard_matrix <- function(celltype_df, colname1, colname2) {
  # make lists of barcodes for each category, named by the category
  id1_list <- split(celltype_df$barcodes, celltype_df[[colname1]])
  id2_list <- split(celltype_df$barcodes, celltype_df[[colname2]])
  
  # create the grid of comparisons
  cross_df <- tidyr::expand_grid(id1 = names(id1_list), id2 = names(id2_list))
  
  # calculate a single Jaccard index for each combination using split lists & ids
  jaccard_scores <- cross_df |>
    purrr::pmap_dbl(\(id1, id2){
      jaccard(id1_list[[id1]], id2_list[[id2]])
    })
  
  # add scores to the comparison grid and convert to matrix
  jaccard_matrix <- cross_df |>
    dplyr::mutate(jaccard = jaccard_scores) |>
    # convert to matrix
    tidyr::pivot_wider(
      names_from = "id2",
      values_from = "jaccard"
    ) |>
    tibble::column_to_rownames(var = "id1") |>
    as.matrix()
  
  return(jaccard_matrix)
}
