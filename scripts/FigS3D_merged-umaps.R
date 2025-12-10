# This script is used to generate a small version of the merged UMAP for Figure S3D
# By default, it used 4 threads.

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)
library(optparse)


option_list <- list(
  make_option(
    "--threads",
    default = 4,
    type = "integer",
    help = "Number of threads to use"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))


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
      panel.background = element_rect(colour = "black", linewidth = 0.5),
      aspect.ratio = 1
    )
)

# Set up -----------------------------------------------------------------------

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

libraries_to_plot <- c(
  "SCPCL000050",
  "SCPCL000697",
  "SCPCL000698",
  "SCPCL000705"
)

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCP000003")

# define file paths to sce files
rds_files <- list.files(
  path = local_results_dir,
  pattern = "_processed.rds",
  full.names = TRUE
) |>
  purrr::set_names(\(x) {
    stringr::str_split_i(basename(x), "_", 1)
  })

# Ensure we have the correct libraries
rds_files <- rds_files[libraries_to_plot]
stopifnot(
  "Libraries needed to make figure S3D are not all present" = length(setdiff(
    libraries_to_plot,
    names(rds_files)
  )) ==
    0
)

# define output file paths
png_dir <- here::here("figures", "pngs")
output_png_file <- file.path(png_dir, "FigS3D_merged-umaps.png")

pdf_dir <- here::here("figures", "pdfs")
output_pdf_file <- file.path(pdf_dir, "FigS3D_merged-umaps.pdf")


# Read, merge, and process objects ------------------------------------------------

# Read in files and merge into a single object
sce_list <- rds_files |>
  purrr::map(
    \(rds_file) {
      sce <- readr::read_rds(rds_file)

      # Add library id to coldata before merging
      sce$library_id <- metadata(sce)$library_id

      return(sce)
    }
  )
# note that this will produce warnings from rowData merging, but these are ok; this is not relevant to the plot
merged_sce <- do.call(combineCols, unname(sce_list))

# clean up for space
rm(sce_list)

# Perform batch-aware PCA and calculate UMAP
# Code adapted from `scpca-nf`: https://github.com/AlexsLemonade/scpca-nf/blob/5cbacbf5be11e8b56c6595bc434eb4e85a483865/bin/merge_sces.R#L321

batch_column <- merged_sce$library_id

# model gene variance & identify subset of variable genes
hvg_list <- scran::modelGeneVar(
  merged_sce,
  block = batch_column,
  BPPARAM = bp_param
) |>
  scran::getTopHVGs(n = 2000)

# perform multi batch PCA on merged object
multi_pca <- batchelor::multiBatchPCA(
  merged_sce,
  subset.row = hvg_list,
  batch = batch_column,
  preserve.single = TRUE,
  BPPARAM = bp_param
)

# add PCA results to merged SCE object
reducedDim(merged_sce, "PCA") <- multi_pca[[1]]

# add UMAP
merged_sce <- scater::runUMAP(
  merged_sce,
  dimred = "PCA"
)

# Plot the UMAP -------------
umap_df <- data.frame(
  library_id = batch_column,
  UMAP1 = reducedDim(merged_sce, "UMAP")[, 1],
  UMAP2 = reducedDim(merged_sce, "UMAP")[, 2]
)

umap_plot <- ggplot(
  umap_df,
  aes(x = UMAP1, y = UMAP2, color = library_id)
) +
  # set points for all "other" points
  geom_point(
    data = dplyr::select(
      umap_df,
      -c("library_id")
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
