# This script is used to generate a small version of the merged UMAP for Figure 3

# load project
renv::load()

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

# Set default ggplot theme
theme_set(
  theme_classic() +
    theme(
      #plot.margin = margin(rep(20, 4)),
      strip.background = element_rect(fill = "transparent", linewidth = 0.5),
      # no axis ticks or labels
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # add a square around each of the plots
      panel.background = element_rect(colour = "black", linewidth=0.5),
      aspect.ratio = 1
    )
)

# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCP000003")

# define file path to merged object 
merged_sce_file <- file.path(local_results_dir, "SCPCP000003_merged.rds")

# read in merged sce objects 
merged_sce <- readr::read_rds(merged_sce_file)

# define output file paths 
png_dir <- here::here("figures", "pngs") 
output_png_file <- file.path(png_dir, "Fig3D_merged-umaps.png")

pdf_dir <- here::here("figures", "pdfs") 
output_pdf_file <- file.path(pdf_dir, "Fig3D_merged-umaps.pdf")


# Prep for plotting ------------------------------------------------------------

libraries_to_plot <- c("SCPCL000050",
                       "SCPCL000697",
                       "SCPCL000698",
                       "SCPCL000705")

# extract coldata as data frame to use to create tables and UMAPs
coldata_df <- scuttle::makePerCellDF(merged_sce, use.dimred = "UMAP") |>
  dplyr::rename(
    UMAP1 = UMAP.1,
    UMAP2 = UMAP.2
  ) |> 
  dplyr::filter(library_id %in% libraries_to_plot)

# UMAPs ------------------------------------------------------------------------

umap_plot <- ggplot(
  coldata_df,
  aes(x = UMAP1, y = UMAP2, color = library_id)
) +
  # set points for all "other" points
  geom_point(
    data = dplyr::select(
      coldata_df, -c("library_id")
    ),
    color = "gray80",
    alpha = 0.5,
    size = 0.1
  ) +
  # set points for desired cell type
  geom_point(
    alpha = 0.5,
    size = 0.1,
    color = "firebrick3"
  ) +
  facet_wrap(
    vars(library_id), 
    ncol = 4 # make one line of UMAPs
  )

# export files 
ggsave(output_png_file, umap_plot, width = 10, height = 7)
ggsave(output_pdf_file, umap_plot, width = 10, height = 7)
