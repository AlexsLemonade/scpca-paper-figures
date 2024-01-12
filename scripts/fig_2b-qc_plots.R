# This script is used to generate simplified version of all the plots in the main QC report 

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
      strip.background = element_rect(fill = "transparent"),
      # no background box for legend
      legend.background = element_blank(),
      # no axis ticks or labels
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # add a square around each of the plots
      panel.background = element_rect(colour = "black", linewidth=0.5),
      aspect.ratio = 1
    )
)

set.seed(2024)

# Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "SCPCS000001")

# define file paths to downloaded files 
unfiltered_sce_file <- file.path(local_results_dir, "SCPCL000001_unfiltered.rds")
filtered_sce_file <- file.path(local_results_dir, "SCPCL000001_filtered.rds")
processed_sce_file <- file.path(local_results_dir, "SCPCL000001_processed.rds")

# read in sce objects 
unfiltered_sce <- readr::read_rds(unfiltered_sce_file)
filtered_sce <- readr::read_rds(filtered_sce_file)
processed_sce <- readr::read_rds(processed_sce_file)

# define output file paths 
plots_dir <- here::here("figures", "pngs") 
output_plot_file <- file.path(plots_dir, "Fig2B-mini-qc-plots.png")

# Knee plot --------------------------------------------------------------------

# Create a knee plot using the same code from scpca-nf 
unfiltered_celldata <- data.frame(colData(unfiltered_sce)) |>
  mutate(
    rank = rank(-unfiltered_sce$sum, ties.method = "first"), # using full spec for clarity
    filter_pass = colnames(unfiltered_sce) %in% colnames(filtered_sce)
  ) |>
  select(sum, rank, filter_pass) |>
  filter(sum > 0) # remove zeros for plotting


grouped_celldata <- unfiltered_celldata |>
  mutate(rank_group = floor(rank / 100)) |>
  group_by(rank_group) |>
  summarize(
    med_sum = median(sum),
    med_rank = median(rank),
    pct_passed = sum(filter_pass) / n() * 100
  )

top_celldata <- unfiltered_celldata |>
  filter(rank <= 50) |>
  mutate(filter_pct = ifelse(filter_pass, 100, 0))

knee_plot <- ggplot(grouped_celldata, aes(x = med_rank, y = med_sum, color = pct_passed)) +
  geom_point(
    mapping = aes(x = rank, y = sum, color = filter_pct),
    data = top_celldata,
    alpha = 0.5,
    size = 0.5
  ) +
  geom_line(linewidth = 1, lineend = "round", linejoin = "round") +
  scale_x_log10(labels = scales::label_number(accuracy = 1)) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  scale_color_gradient2(
    low = "grey70",
    mid = "forestgreen",
    high = "darkgreen",
    midpoint = 50
  ) +
  labs(
    x = "Rank",
    y = "Total UMI count",
    color = ""
  ) +
  theme(
    legend.position = c(.95, .9),
    aspect.ratio = 1
    
  ) +
  # remove legend labels
  guides(color = guide_colorbar(label = FALSE, barwidth = 0.5, barheight = 2))

# Cell read metrics ------------------------------------------------------------

# create metrics of reads/ genes detected from qc report 
filtered_celldata <- data.frame(colData(filtered_sce))

cell_metrics_plot <- ggplot(
  filtered_celldata,
  aes(
    x = sum,
    y = detected,
    color = subsets_mito_percent
  )
) +
  geom_point(alpha = 0.3, size = 0.2) +
  scale_color_viridis_c(limits = c(0, 100)) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  labs(
    x = "Total UMI count",
    y = "Number of genes detected",
    color = ""
  ) +
  theme(
    legend.position = c(.95, .9)
    
  ) +
  guides(color = guide_colorbar(label = FALSE, barwidth = 0.5, barheight = 2))

# miQC metrics -----------------------------------------------------------------

# create miQC plot 
filtered_sce$prob_compromised <- NULL
miQC_model <- metadata(filtered_sce)$miQC_model

miQC_plot <- miQC::plotModel(filtered_sce, model = miQC_model) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
  labs(
    x = "Number of genes detected",
    y = "Percent reads mitochondrial",
    color = ""
  ) +
  theme(
    legend.position = c(.95, .9)
    
  ) +
  guides(color = guide_colorbar(label = FALSE, barwidth = 0.5, barheight = 2))

# set line thickness
line_aes <- list(linewidth = 0.5, alpha = 0.8)
miQC_plot$layers[[2]]$aes_params <- line_aes
miQC_plot$layers[[3]]$aes_params <- line_aes

# set point size and transparency
miQC_plot$layers[[1]]$aes_params <- list(size=0.2, alpha = 0.5)


# Filter low quality -----------------------------------------------------------

# grab cutoffs and filtering method from processed sce
processed_meta <- metadata(processed_sce)
min_gene_cutoff <- processed_meta$min_gene_cutoff
filter_method <- processed_meta$scpca_filter_method

# add column to coldata labeling cells to keep/remove based on filtering method
filtered_coldata_df <- colData(filtered_sce) |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode")

# plot showing low quality cells that are filtered out
filtered_plot <- ggplot(filtered_coldata_df, aes(x = detected, y = subsets_mito_percent, color = scpca_filter)) +
  geom_point(alpha = 0.5, size = 0.2) +
  labs(
    x = "Number of genes detected",
    y = "Mitochondrial percentage",
    color = ""
  ) +
  theme(
    legend.position = c(.8, .9),
    legend.key = element_blank()
    
  ) + 
  guides(color = guide_legend(override.aes = list(size = 2)))

# UMAP -------------------------------------------------------------------------

# umap showing total number of genes detected
umap_plot <- scater::plotUMAP(
  processed_sce,
  point_size = 0.3,
  point_alpha = 0.5,
  colour_by = "detected"
) +
  scale_color_viridis_c() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(colour = "black", linewidth = 0.5),
        aspect.ratio = 1)

# HVGs -------------------------------------------------------------------------

# select top genes to plot
all_hvg <- processed_meta$highly_variable_genes
# select 4 of the top 5 genes that have short gene names 
# skip the first one since it doesn't have a gene symbol 
top_genes <- all_hvg[2:5]

# grab expression for top genes from counts
var_gene_exp <- logcounts(processed_sce[top_genes, ]) |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode")

# grab rowdata as data frame to later combine with gene expression data
# rowdata contains the mapped gene symbol so can use that to label plots instead of ensembl gene id from rownames
rowdata_df <- rowData(processed_sce) |>
  as.data.frame() |>
  tibble::rownames_to_column("ensembl_id") |>
  select(ensembl_id, gene_symbol) |>
  filter(ensembl_id %in% top_genes) |>
  mutate(
    gene_symbol = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id),
    ensembl_id = factor(ensembl_id, levels = top_genes)
  ) |>
  arrange(ensembl_id) |>
  mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol))


# extract umap embeddings as a dataframe to join with gene expression and coldata for plotting
umap_df <- reducedDim(processed_sce, "UMAP") |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode")

# combine gene expression with coldata, umap embeddings, and rowdata and create data frame to use for plotting
coldata_df <- colData(processed_sce) |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode") |>
  # combine with gene expression
  left_join(var_gene_exp, by = "barcode") |>
  # combine with umap embeddings
  left_join(umap_df, by = "barcode") |>
  # combine all genes into a single column for easy faceting
  tidyr::pivot_longer(
    cols = starts_with("ENSG"),
    names_to = "ensembl_id",
    values_to = "gene_expression"
  ) |>
  # join with row data to add in gene symbols
  left_join(rowdata_df)

# faceted plot for hvgs
hvg_plot <- ggplot(coldata_df, aes(x = UMAP1, y = UMAP2, color = gene_expression)) +
  geom_point(alpha = 0.1, size = 0.001) +
  facet_wrap(vars(gene_symbol)) +
  scale_color_viridis_c() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 8)
  ) 

# Combine plots ----------------------------------------------------------------

plot_list <- list(knee_plot,
                  cell_metrics_plot,
                  miQC_plot,
                  filtered_plot,
                  umap_plot, 
                  hvg_plot)
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 3) &
  theme(text = element_text(size = 10))

# save files 
ggsave(output_plot_file, plot = combined_plot, width = 8.5, height = 5.5, units = "in")
