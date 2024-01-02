# Generates a bar plot summarizing the diagnoses type in the ScPCA Portal 
# The resulting plot is faceted by the diagnoses group 
# Individual bars have different patterns based on the disease timing 

# load any libaries 
library(ggplot2)

# Define functions -------------------------------------------------------------

#' Create individual stacked bar plot panels with diagnosis on y-axis and patterns showing disease timing
#'
#' @param diagnosis_plot_df DataFrame that contains diagnosis, diagnosis_group, and standardized_disease_timing
#'   as columns 
#' @param diagnosis_group Which group of diagnosis to include in the plot, must be a value in diagnosis_group column
#' @param group_color What color to use for the bar plot
#'
#' @return ggplot with diagnosis on y-axis, number of samples on x-axis, and bars filled in with patterns to 
#'   designate disease timing group 

barplot_panel <- function(diagnosis_plot_df, diagnosis_group, group_color){
  
  # pull out the group to be included in this panel
  group_to_plot <- diagnosis_plot_df |> 
    dplyr::filter(diagnosis_group %in% {{diagnosis_group}}) 
  
  # get a count for total diagnosis
  # use this later to sort the diagnosis
  diagnosis_count <- group_to_plot |> 
    dplyr::count(diagnosis) 
  
  # create dataframe with diagnosis, diagnosis count, and disease timing for plotting
  plot_df <- group_to_plot |> 
    dplyr::left_join(diagnosis_count) |>
    dplyr::mutate(diagnosis = forcats::fct_reorder(diagnosis, n))
  
  
  diagnosis_plot <- ggplot(plot_df, aes(y = diagnosis)) + 
    # set pattern aesthetics 
    ggpattern::geom_bar_pattern(aes(pattern = standardized_disease_timing),
                                pattern_density = 0.04,
                                fill = group_color,
                                color = "black",
                                pattern_color = "black",
                                pattern_spacing = 0.02,
                                pattern_orientation = 'vertical') +
    theme_classic() + 
    theme(text = element_text(size = 10),
          # make sure legend is big enough to see the patterns
          legend.key.size = unit(1, 'cm'),
          # don't print out the legend for the individual panels
          legend.position = 'none') +
    labs(
      x = "",
      y = "",
      title = diagnosis_group,
      color = "Disease timing"
    ) +
    # create a custom pattern scale so that we keep the patterns the same for each disease timing
    ggpattern::scale_pattern_manual(values = c(
      "Initial diagnosis" = 'circle', 
      "Post-mortem" = 'stripe', 
      "Progressive" = 'crosshatch', 
      "Recurrence" = 'none',
      "Unknown" = 'wave'), 
      drop = FALSE) 
  
  return(diagnosis_plot)
}


# Set up -----------------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# path to sample metadata stored on S3
sample_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-sample-metadata.tsv'

# create a directory to store S3 files 
local_s3_dir <- file.path(root_dir, "s3_files")
fs::dir_create(local_s3_dir)

# all metadata files 
sample_info_dir <- file.path(root_dir, "sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
diagnosis_groupings_file <- file.path(sample_info_dir, "diagnosis-groupings.tsv")
disease_timing_file <- file.path(sample_info_dir, "disease-timing.tsv")

# output files 
plots_dir <- file.path(root_dir, "figures", "pngs")
output_plot_file <- file.path(plots_dir, "Fig1-sample-summary.png")

# Prep sample metadata ------------------------------------------------------

# read in project whitelist and grouping metadata 
project_whitelist <- readLines(project_whitelist_file)
diagnosis_groupings_df <- readr::read_tsv(diagnosis_groupings_file) |>
  dplyr::select(submitted_diagnosis, diagnosis_group)
disease_timing_df <- readr::read_tsv(disease_timing_file)

# Copy sample metadata 
local_sample_metadata <- file.path(local_s3_dir, "scpca-sample-metadata.tsv")
sync_call <- paste('op run -- aws s3 cp', sample_metadata_s3, local_sample_metadata, '--recursive', sep = " ")
system(sync_call)

# read in sample metadata and filter to only projects in whitelist
sample_metadata_df <- readr::read_tsv(local_sample_metadata) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist)

# Join sample metadata with diagnosis grouping and disease timing 
diagnosis_plot_df <- sample_metadata_df |> 
  dplyr::select(scpca_sample_id, diagnosis, disease_timing) |> 
  dplyr::left_join(diagnosis_groupings_df, by = c("diagnosis" = "submitted_diagnosis")) |> 
  dplyr::left_join(disease_timing_df, by = c("disease_timing" = "submitted_disease_timing")) |> 
  # right now let's remove this group
  dplyr::filter(diagnosis_group != "Non-cancerous") |> 
  # make it a factor so we can use ggpattern 
  dplyr::mutate(standardized_disease_timing = as.factor(standardized_disease_timing))


# get list of all groups 
diagnosis_groups <- diagnosis_plot_df |> 
  dplyr::pull(diagnosis_group) |>
  unique() |>
  # set desired order 
  forcats::fct_relevel("Brain and CNS",
                       "Leukemia", 
                       "Sarcoma",
                       "Other solid tumors")
# colors 
# TODO: move these into a separate tsv file to be read in 
#names(diagnosis_groups) <- c("#201158", "#005E5E", "#578B21", "#E89E6B")
names(diagnosis_groups) <- c("#DF536B", "#61D04F", "#2297E6", "#28E2E5")

# Plot -------------------------------------------------------------------------

# Create list of plots that will be combined into one figure 
all_panels <- diagnosis_groups |> 
  purrr::imap(\(group, group_color){
    barplot_panel(diagnosis_plot_df = diagnosis_plot_df,
                  diagnosis_group = group, 
                  group_color = group_color)
  })

# grab one legend from list of plots to use as the common legend
legend_to_use <- ggpubr::get_legend(all_panels[[1]])

# combine into one plot using a common legend 
combined_plot <- ggpubr::ggarrange(plotlist = all_panels,
                                   common.legend = TRUE, 
                                   legend = "top",
                                   legend.grob = legend_to_use, 
                                   ncol = 2, 
                                   nrow = 2, 
                                   labels = "AUTO",
                                   align = 'hv',
                                   hjust = -5) +
  # make sure white background encompasses all panels and legend 
  ggpubr::bgcolor("white")

# add x and y axis labels 
combined_plot <- ggpubr::annotate_figure(combined_plot, 
                                         bottom = "Number of samples",
                                         left = "Diagnosis")
  

# save plot 
ggsave(output_plot_file, width = 15)
