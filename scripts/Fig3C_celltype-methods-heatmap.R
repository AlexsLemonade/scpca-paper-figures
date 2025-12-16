# This script is generates a heatmap comparing automated method cell type annotations to consensus cell type annotations

# load project
renv::load()

library(SingleCellExperiment)

# Set heatmap padding option
ComplexHeatmap::ht_opt(TITLE_PADDING = grid::unit(0.2, "in"))

# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCS000001")

# define file paths to downloaded files 
processed_sce_file <- file.path(local_results_dir, "SCPCL000001_processed.rds")

# read in sce objects 
processed_sce <- readr::read_rds(processed_sce_file)

# define output file paths 
figure_dir <- here::here("figures", "pdfs")
output_file <- file.path(figure_dir, "Fig3C_celltype-comparison-heatmap.pdf")

# source in helper functions for plotting
function_file <- here::here("scripts", "utils", "celltype-plot-helper-functions.R")
source(function_file)

# Prep data for plotting -------------------------------------------------------
top_n_consensus <- 7 # show the top 7 consensus cell types
n_cells <- 3 # we'll only show cell types with >= 3 cells
unknown_labels <- c("Unclassified cell", "Unknown", "other") # collapse into one unknown grouping

# data frame of cell types
celltype_df <- colData(processed_sce) |>
  as.data.frame() |>
  dplyr::select(
    barcodes,
    singler = singler_celltype_annotation, 
    cellassign = cellassign_celltype_annotation, 
    scimilarity = scimilarity_celltype_annotation, 
    consensus = consensus_celltype_annotation
  ) |>
  tibble::as_tibble()

# identify top n consensus to show
top_consensus <- celltype_df |>
  dplyr::count(consensus)|>
  dplyr::slice_max(n, n = top_n_consensus) |>
  dplyr::pull(consensus)

# prepare data frame to make heatmap from
heatmap_celltype_df <- celltype_df |>
  dplyr::filter(consensus %in% top_consensus) |>
  # pivot for filtering down and recoding celltypes
  tidyr::pivot_longer(
    -barcodes, 
    names_to = "method", 
    values_to = "celltype"
  ) |>
  # only plot cell types with at least n_cells cells
  dplyr::add_count(method, celltype) |>
  dplyr::filter(n >= n_cells) |>
  dplyr::select(-n) |>
  dplyr::mutate(
    # use the same label for all unknowns 
    celltype = ifelse(
      is.na(celltype) | celltype %in% unknown_labels, "unknown", celltype
    ), 
    # lowercase for manuscript consistency
    celltype = stringr::str_to_lower(celltype)
  ) |> 
  # back to wide for making the matrix
  tidyr::pivot_wider(
    names_from = "method", 
    values_from = "celltype"
  ) |>
  # set factor levels for each method separately, in order of frequency with "unknown" at the end
  dplyr::mutate(
    dplyr::across(
      -barcodes,
      \(x) forcats::fct_infreq(x) |> forcats::fct_relevel("unknown", after = Inf)
    )
  ) 

# Update consensus labels only to include counts
consensus_labels <- heatmap_celltype_df |>
  dplyr::count(consensus) |>
  dplyr::arrange(consensus) |>
  dplyr::mutate( 
    consensus_n = glue::glue("{consensus} (n = {n})")
  ) |>
  dplyr::pull(consensus_n)
heatmap_celltype_df$consensus <- factor(
  heatmap_celltype_df$consensus, 
  levels = levels(heatmap_celltype_df$consensus), 
  labels = consensus_labels
) 


# prepare list of jaccard matrices to use in heatmap
jaccard_list <- method_names |>
  purrr::map(\(name) {
    make_jaccard_matrix(
      heatmap_celltype_df,
      "consensus",
      name
    )
  })

# Make and export the heatmap -------------------------


################## TODO: PICK A VERSION TO PROCEED WITH ################

# show methods in alphabetical order
method_names <- c(
  "CellAssign annotations" = "cellassign", 
  "SCimilarity annotations" = "scimilarity",
  "SingleR annotations" = "singler"
)

## VERTICAL VERSION
heatmap_list <- jaccard_list |>
  purrr::imap(
    \(mat, name) {
      ComplexHeatmap::Heatmap(
        t(mat), 
        col = circlize::colorRamp2(c(0, 1), colors = c("white", "darkslateblue")),
        border = TRUE, # each heatmap gets its own outline
        ## Row parameters
        cluster_rows = FALSE,
        row_title_gp = grid::gpar(fontsize = 10),
        row_title_side = "left",
        row_title = name,
        row_names_side = "right",
        row_names_gp = grid::gpar(fontsize = 10),
        ## Column parameters
        cluster_columns = FALSE,
        column_title = "Consensus cell type annotations",
        column_title_gp = grid::gpar(fontsize = 10),
        column_names_side = "bottom",
        column_names_gp = grid::gpar(fontsize = 10),
        ## Legend parameters
        heatmap_legend_param = list(
          title = "Jaccard index",
          direction = "horizontal",
          legend_width = grid::unit(1.5, "in")
        ),
        # make sure we only have 1 legend
        show_heatmap_legend = name == "SingleR annotations",
      )
    }
  ) |>
  # concatenate TBD into HeatmapList object
  purrr::reduce(ComplexHeatmap::`%v%`) # use + for horizontal

final_heatmap <- ComplexHeatmap::draw(
  heatmap_list_vertical,
  heatmap_legend_side = "bottom"
)


## HORIZONTAL VERSION
heatmap_list <- jaccard_list |>
  purrr::imap(
    \(mat, name) {
      ComplexHeatmap::Heatmap(
        mat, 
        col = circlize::colorRamp2(c(0, 1), colors = c("white", "darkslateblue")),
        border = TRUE, # each heatmap gets its own outline
        ## Row parameters
        cluster_rows = FALSE,
        row_title_gp = grid::gpar(fontsize = 10),
        row_title_side = "left",
        row_title = "Consensus cell type annotations",
        row_names_side = "right",
        row_names_gp = grid::gpar(fontsize = 10),
        ## Column parameters
        cluster_columns = FALSE,
        column_title = name,
        column_title_gp = grid::gpar(fontsize = 10),
        column_names_side = "bottom",
        column_names_gp = grid::gpar(fontsize = 10),
        ## Legend parameters
        heatmap_legend_param = list(
          title = "Jaccard index",
          direction = "horizontal",
          legend_width = grid::unit(1.5, "in")
        ),
        # make sure we only have 1 legend
        show_heatmap_legend = name == "SingleR annotations",
      )
    }
  ) |>
  # concatenate TBD into HeatmapList object
  purrr::reduce(`+`) 
###############################################################################


# save heatmap to pdf
pdf(output_file, width = 9, height = 9, useDingbats = FALSE) # TODO: determine dimensions
ComplexHeatmap::draw(
  heatmap_list,
  heatmap_legend_side = "bottom"
)
dev.off()
