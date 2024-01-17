# This is the script used to make the benchmarking figures that are also 
# found in scpca-docs comparing Alevin-fry to Cell Ranger 
# Prior to running this script, you will need to download the input data using 
# `figure_setup/sync-data-files.R` and the reference files using `figure_setup/sync-reference-files.R`

renv::load()

library(ggplot2)
library(SingleCellExperiment)

# Set up -----------------------------------------------------------------------

# list of libraries 
single_cell <- c("SCPCR000003", "SCPCR000126", "SCPCR000127")
single_nuclei <- c("SCPCR000220", "SCPCR000221", "SCPCR000495")
all_libs <- c(single_cell, single_nuclei)

# directory with benchmarking results
benchmark_dir <- here::here("s3_files", "benchmarking_results")

# list of alevin fry directories to read in as SCE
af_dir <- file.path(benchmark_dir, "alevin-fry")
af_dir_list <- file.path(af_dir, all_libs) |> 
  purrr::set_names(all_libs)

# list of cellranger directories to read in as SCE 
cellranger_dir <- file.path(benchmark_dir, "cellranger")
cellranger_dir_list <- file.path(cellranger_dir, all_libs) |> 
  purrr::set_names(all_libs)

# read in mito genes 
mito_file <- here::here("s3_files", "reference_files", "Homo_sapiens.GRCh38.104.mitogenes.txt")
mito_genes <- readLines(mito_file)

# Create SCE objects -----------------------------------------------------------

# read in alevin-fry output as SCEs
af_sces <- af_dir_list |> 
  purrr::map(
    \(quant_dir){
      sce <- scpcaTools::read_alevin(
        quant_dir = quant_dir,
        fry_mode = TRUE) |> 
        # filter empty droplets 
        scpcaTools::filter_counts()
    }
  )

# read in cell ranger output as SCEs 
cellranger_sces <- cellranger_dir_list |> 
  purrr::map(scpcaTools::read_cellranger) 

# add per cell and per feature info to each sce 
all_sces <- list(af_sces, cellranger_sces) |> 
  purrr::set_names(c("Alevin-fry", "Cell Ranger")) |> 
  purrr::map(\(list){
    list |> 
      purrr::map(\(sce) {
        sce |> 
          scuttle::addPerCellQCMetrics(subsets = list(mito = mito_genes[mito_genes %in% rownames(sce)])) |>
          scuttle::addPerFeatureQCMetrics()
      })
  })

# Prep for plotting ------------------------------------------------------------

# create a data frame with coldata info for each tool, run id combo 
coldata_df <- all_sces |> 
  purrr::map_df(~ purrr::map_df(.x, scpcaTools::coldata_to_df, .id = "run_id"),
    .id = "tool"
    ) |> 
  # create new columns with seq unit 
  dplyr::mutate(seq_unit = dplyr::case_when(run_id %in% single_cell ~ "Cell",
                                            run_id %in% single_nuclei ~ "Nuclei"),
                # create title for plots 
                plot_id = glue::glue("{run_id}-{seq_unit}"))

# filter for cells that are found in both af + cellranger
cell_counts <- coldata_df |>  
  dplyr::count(cell_id, run_id)

common_cells <- cell_counts |>
  dplyr::filter(n == 2) |>
  dplyr::pull(cell_id)

coldata_common <- coldata_df |>
  dplyr::filter(cell_id %in% common_cells) 

# Plot cell metrics ------------------------------------------------------------

# create combined UMI per cell plot 
umi_plot <- ggplot(coldata_common, aes(x = sum, color = tool)) + 
  geom_density() + 
  facet_wrap(vars(plot_id)) +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "top") +
  labs(x = expression(paste(Log[10], " total UMI per cell")),
       y = "") +
  #TODO: Replace with palette colors 
  scale_color_manual(values = c("#D55E00", "#0072B2"))

genes_plot <- ggplot(coldata_common, aes(x = detected, color = tool)) + 
  geom_density() + 
  facet_wrap(~ plot_id) +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "top") +
  labs(x = expression(paste(Log[10], " total genes detected per cell")),
       y = "") +
  scale_color_manual(values = c("#D55E00", "#0072B2"))

cell_metrics_plots <- patchwork::wrap_plots(list(umi_plot, genes_plot), ncol = 1, guides = 'collect')

# prep row data ----------------------------------------------------------------

## Mean gene expression correlation  
# grab rowdata from filtered sces
rowdata_df <- all_sces |> 
  purrr::map_df(
    ~ purrr::map_df(.x, scpcaTools::rowdata_to_df, .id = "run_id"),
    .id = "tool"
  ) |>
  dplyr::mutate(seq_unit = dplyr::case_when(run_id %in% single_cell ~ "Cell",
                                            run_id %in% single_nuclei ~ "Nuclei"),
                plot_id = glue::glue("{run_id}-{seq_unit}"))

# remove genes with low detection 
rowdata_cor <- rowdata_df |>
  dplyr::filter(detected >= 5.0) |> 
  # spread table to put mean expression for Alevin-fry and Cell ranger in its own columns for plotting
  dplyr::select(tool, gene_id, plot_id, mean) |>
  # spread the mean expression stats to one column per caller
  tidyr::pivot_wider(id_cols = c(gene_id, plot_id),
                     names_from = c("tool"),
                     values_from = mean) |>
  # drop rows with NA values to ease correlation calculations
  tidyr::drop_na()

# Plot gene metrics ------------------------------------------------------------

# plot correlation of alevin-fry to Cellranger for each sample 
gene_exp_plot <- ggplot(rowdata_cor, aes(x = `Alevin-fry`, y = `Cell Ranger`)) +
  geom_point(size = 0.5, alpha = 0.1) + 
  facet_wrap(~ plot_id) + 
  scale_x_log10(labels = scales::label_number()) + 
  scale_y_log10(labels = scales::label_number()) + 
  labs(x = "Cell Ranger mean gene expression", y = "Alevin-fry mean gene expression") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 6),
        strip.background = element_rect()) +
  ggpubr::stat_cor(aes(label = after_stat(r.label)), method = "pearson", size = 2)

cell_metrics_plots / gene_exp_plot
