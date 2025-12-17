# This script generates the panels of S4 which shows detailed cell type annotations 
# and CNV information for a single sample, SCPCS000049/SCPCL000049

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

# set default ggplot theme for UMAPs
theme_set(
  theme_classic() +
    theme(
      strip.background = element_rect(fill = "transparent", linewidth = 0.5),
      # no axis ticks or labels
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # add a square around each of the plots
      panel.background = element_rect(colour = "black", linewidth = 0.5),
      aspect.ratio = 1,
      # remove boxes around legends 
      legend.key = element_blank()
    )
)
# Set up -----------------------------------------------------------------------

sce_file <- here::here("s3_files", "SCPCS000049", "SCPCL000049_processed.rds")

# define output files
tcell_umap_file <- here::here("figures", "pdfs", "FigS4A_tcell-umaps.pdf")
consensus_cnv_umap_file <- here::here("figures", "pdfs", "FigS4B_consensus-cnv-umaps.pdf")

# get cell types and umap embeddings for plotting
umap_df <- readr::read_rds(sce_file) |>
  # convert to data frame immediately for space
  scuttle::makePerCellDF(
    use.coldata = c(
      "barcodes",
      "consensus_celltype_annotation",
      "singler_celltype_annotation",
      "cellassign_celltype_annotation",
      "scimilarity_celltype_annotation",
      "infercnv_total_cnv"
    ),
    use.dimred = "UMAP"
  ) |>
  # remove rownames
  tibble::as_tibble()

# Faceted UMAP - panel A -------------------------------------------------------

# list of cell type columns to lump and plot
celltype_columns <- c(
  "Consensus" = "consensus_celltype_annotation",
  "SingleR" = "singler_celltype_annotation",
  "CellAssign" = "cellassign_celltype_annotation",
  "SCimilarity"= "scimilarity_celltype_annotation"
)

# make a separate faceted plot for each cell type method
# show the top 3 T cell types for each method and then lump all other T cells
all_plot_list <- celltype_columns |> 
  purrr::imap(\(celltype_column, method) {
    
    # lump the T cell groups
    celltype_df <- umap_df |> 
      dplyr::mutate(
        t_cell_group = dplyr::if_else(
          # keep T cells with the original annotation and leave all others as NA
          stringr::str_detect(!!sym(celltype_column), "T[- ]"),
          !!sym(celltype_column),
          NA_character_
        ) |> 
          # get just the top 3 T cell types and then lump all others together
          forcats::fct_lump_n(n = 3, other_level = "Other T cells") |> 
          stringr::str_wrap(25) |> # wrap the long names
          forcats::fct_infreq() |> 
          forcats::fct_relevel("Other T cells", after = Inf)
      )
    
    # update the labels to include counts, drop NA since that's not going to be a level
    celltype_labels <- celltype_df |> 
      dplyr::count(t_cell_group) |> 
      dplyr::mutate(
        facet_label = glue::glue("{t_cell_group} (n = {n})")
      ) |> 
      tidyr::drop_na() |> 
      dplyr::pull(facet_label)
    
    # add a new column for labels and filter
    plot_df <- celltype_df |> 
      dplyr::mutate(
        facet_label = factor(t_cell_group, levels = levels(t_cell_group), labels = celltype_labels)
      ) |> 
      # filter to only keep panels that have T cells
      # otherwise we get an NA panel 
      dplyr::filter(!is.na(t_cell_group))
    
    ggplot(
      plot_df,
      aes(x = UMAP.1, y = UMAP.2, color = t_cell_group)
    ) +
      # set points for all "other" points
      geom_point(
        data = dplyr::select(
          celltype_df, -t_cell_group
        ),
        color = "gray80",
        alpha = 0.5,
        size = 0.3
      ) +
      # set points for desired cell type
      geom_point(size = 0.3, alpha = 0.5, color = "firebrick3") +
      # remove axis numbers and background grid
      scale_x_continuous(labels = NULL, breaks = NULL) +
      scale_y_continuous(labels = NULL, breaks = NULL) +
      labs(
        x = NULL,
        y = method # label with celltype method
      ) +
      # facet by t cell type
      facet_grid(cols = vars(facet_label)) +
      theme(
        # Ensure plot margins are tight so the axis labels can be as close as possible
        plot.margin = margin(1, 0, 0, 1)
      )
    
  }) |> 
  # combine so that each method gets its own row
  patchwork::wrap_plots(ncol = 1)

# annotate the plot to have a single combined x and y axes for UMAP1 and UMAP2
y_label_global <- wrap_elements(grid::textGrob("UMAP2", rot = -90))
x_label_global <- wrap_elements(grid::textGrob("UMAP1"))

# wrap main plot otherwise patchwork complains when making the layout
main_grob <- wrap_elements(panel = all_plot_list)

# make the combined plot with x and y axes labels
# use units/null in plot layout to make sure the main plot takes up all the space and x axis label is small
combined_plot <- (main_grob + y_label_global + plot_layout(widths = unit(c(1.5, 1), c("null", "lines")))) / 
  x_label_global + 
  plot_layout(heights = unit(c(5, 1), c("null", "lines")))

ggsave(tcell_umap_file, combined_plot, width = 10, height = 10)

# Consensus/inferCNV UMAPs -----------------------------------------------------

# lump consensus cell types 
umap_df <- umap_df |> 
  dplyr::mutate(
    consensus_lumped = forcats::fct_lump_n(consensus_celltype_annotation, n = 7, other_level = "All remaining cell types") |> 
      #stringr::str_wrap(25) |> # wrap the long names
      forcats::fct_infreq() |> 
      forcats::fct_relevel("Unknown", "All remaining cell types", after = Inf)
  )

# get total cell number labels 
celltype_labels <- umap_df |> 
  dplyr::count(consensus_lumped) |> 
  dplyr::mutate(
    celltype_label = glue::glue("{consensus_lumped} (n = {n})")
  ) |> 
  dplyr::pull(celltype_label)

# add a new column for labels
umap_df <- umap_df |> 
  dplyr::mutate(
    celltype_label = factor(consensus_lumped, levels = levels(consensus_lumped), labels = celltype_labels)
  )

# UMAP showing consensus cell types
celltype_plot <- ggplot(umap_df, aes(x = UMAP.1, y = UMAP.2, color = celltype_label)) +
  geom_point(
    alpha = 0.5,
    size = 0.3
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Consensus cell type annotation"
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 1.5))
  )

# cnv umap
cnv_plot <- ggplot(umap_df) +
  aes(x = UMAP.1, y = UMAP.2, color = infercnv_total_cnv) +
  geom_point(alpha = 0.25, size = 0.25) +
  scale_color_viridis_c(na.value = "grey70") +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Total CNV"
  ) 

# combine and export
combined_plot <- celltype_plot + cnv_plot

ggsave(consensus_cnv_umap_file, combined_plot, width = 10, height = 8)
