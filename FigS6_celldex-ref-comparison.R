# This script is used to generate a supplemental figure comparing celldex refs across samples 

# load project
renv::load()

library(SingleCellExperiment)
library(ggplot2)

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      strip.background = element_rect(fill = "transparent"),
      aspect.ratio = 1,
      # no axis ticks or labels
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # add a square around each of the plots
      panel.background = element_rect(colour = "black", linewidth = 0.5),
      # remove boxes around legends 
      legend.key=element_blank()
    )
)

  # Set up -----------------------------------------------------------------------

# folder where any local data lives
local_results_dir <- here::here("s3_files", "celltype_results")

# define file paths to downloaded files 
library_ids <- c(
  "SCPCL000001",
  "SCPCL000002",
  "SCPCL000004",
  "SCPCL000296",
  "SCPCL000297",
  "SCPCL000298",
  "SCPCL000478",
  "SCPCL000484",
  "SCPCL000488"
)

sce_file_paths <- file.path(local_results_dir, glue::glue("{library_ids}_processed.rds")) |> 
  purrr::set_names(library_ids)

# read in all sce objects 
all_sce_objects <- sce_file_paths |>
  purrr::map(readr::read_rds)

# make a list of reference SingleR models and read in
ref_file_list <- c(
  "BlueprintEncodeData" = "BlueprintEncodeData_celldex_1-10-1_model.rds",
  "DatabaseImmuneCellExpressionData" = "DatabaseImmuneCellExpressionData_celldex_1-10-1_model.rds",
  "HumanPrimaryCellAtlasData" = "HumanPrimaryCellAtlasData_celldex_1-10-1_model.rds",
  "MonacoImmuneData" = "MonacoImmuneData_celldex_1-10-1_model.rds"
)
ref_file_dir <- here::here("s3_files", "reference_files", "singler_models")
ref_file_paths <- file.path(ref_file_dir, ref_file_list)
singler_model_list <- ref_file_paths |> 
  purrr::map(readr::read_rds) |> 
  purrr::set_names(names(ref_file_list))

# define path to diagnosis group palette 
diagnosis_group_palette_file <- here::here("palettes", "diagnosis-group-palette.tsv")

# define output file paths 
png_dir <- here::here("figures", "pngs") 
celldex_comparison_png_file <- file.path(png_dir, "FigS6_celldex-ref-comparison.png")

# SingleR ----------------------------------------------------------------------

# for each library, run SingleR for each model 
run_singler <- function(sce, singler_model){
  singler_results <- SingleR::classifySingleR(
    trained = singler_model,
    test = sce,
    fine.tune = TRUE
  )
  return(singler_results)
}

create_delta_median_df <- function(singler_results){
  
  # extract scores into matrix
  singler_scores <- singler_results$scores
  
  # Create data frame for plotting with delta median and the full *non-pruned* cell labels
  delta_median_df <- tibble::tibble(
    delta_median = rowMaxs(singler_scores) - rowMedians(singler_scores),
    # Need to grab the non-pruned label for this plot
    full_labels = singler_results$labels,
    # if pruned.labels are NA ==> low confidence
    # so, negate for this variable:
    confident = !is.na(singler_results$pruned.labels)
  ) |>
    dplyr::mutate(
      confident = ifelse(confident, "High-quality", "Low-quality")
    )
  
  return(delta_median_df)
}

all_delta_df <- all_sce_objects |> 
  # apply to all objects 
  purrr::map(\(sce) {
    # for each singleR model, run SingleR and grab delta median 
    delta_median_df <- singler_model_list |> 
      purrr::map(\(singler_model){
        
        # get singleR results 
        delta_median_df <- run_singler(
          sce = sce,
          singler_model = singler_model
        ) |> 
          # get delta median and confident labels
          create_delta_median_df()
      }) |> 
      # combine by reference 
      dplyr::bind_rows(.id = "celldex_reference")
  }) |>
  # combine by library id 
  dplyr::bind_rows(.id = "library_id")

# Plot -------------------------------------------------------------------------

# define diagnosis groups 
brain_libraries <- c("SCPCL000001", "SCPCL000002", "SCPCL000004")
blood_libraries <- c("SCPCL000296", "SCPCL000297", "SCPCL000298")
sarcoma_libraries <- c("SCPCL000478", "SCPCL000484", "SCPCL000488")

# add diagnosis group
all_delta_df <- all_delta_df |> 
  dplyr::mutate(
    diagnosis_group = dplyr::case_when(
      library_id %in% brain_libraries ~ "Brain and CNS",
      library_id %in% blood_libraries ~ "Leukemia",
      library_id %in% sarcoma_libraries ~ "Sarcoma"
    ),
    plot_title = glue::glue("{library_id}\n{diagnosis_group}")
  ) 

# Subset the data to just confident points for median+/-IQR
delta_median_confident_df <- all_delta_df |>
  dplyr::filter(confident == "High-quality")

# diagnosis group palette 
diagnosis_group_palette <- readr::read_tsv(diagnosis_group_palette_file)

# get list of all colors 
diagnosis_colors <- diagnosis_group_palette$color |> 
  purrr::set_names(diagnosis_group_palette$diagnosis_group)

# define suspension backgrounds to use on facet strips 
backgrounds <- rep(c("Brain and CNS", "Leukemia", "Sarcoma"), each = 3) |>
  purrr::map(\(x) element_rect(fill = diagnosis_colors[[x]]))

# create delta median plot 
all_delta_plot <- ggplot(all_delta_df, aes(x = celldex_reference, y = delta_median, shape = confident, alpha = confident)) +
  ggforce::geom_sina(
    size = 0.5,
    color = "black", # will get applied to all confident points and non-confident outline
    fill = "white", # will apply to non-confident fill only
    position = position_dodge(width = 0.05) # Keep both types of points mostly in line
  ) +
  # Handle points aesthetics:
  #  confident are closed black with alpha = 0.5
  #  not confident are open black with alpha = 1
  scale_shape_manual(values = c(19, 21)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  ggh4x::facet_wrap2(vars(plot_title),
             strip = ggh4x::strip_themed(background_x = backgrounds)) +
  stat_summary(
    data = delta_median_confident_df,
    color = "red",
    geom = "point",
    fun = "median",
    shape = 18,
    size = 2.25,
    alpha = 0.9
  ) +
  labs(
    x = "",
    y = "Delta median statistic",
    shape = "Cell type annotation quality"
  ) +
  guides(
    alpha = FALSE,
    shape = guide_legend(override.aes = list(size = 1.5, alpha = 0.55))
  )  +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)
  )

ggsave(celldex_comparison_png_file, all_delta_plot)
