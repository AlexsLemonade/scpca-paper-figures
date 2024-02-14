# Generates a bar plot summarizing the diagnoses type in the ScPCA Portal 
# The resulting plot is faceted by the diagnoses group 
# Individual bars have different patterns based on the disease timing 

renv::load()

# load any libaries 
library(ggplot2)

# Set up -----------------------------------------------------------------------

# all metadata files 
sample_info_dir <- here::here("sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
diagnosis_groupings_file <- file.path(sample_info_dir, "diagnosis-groupings.tsv")
disease_timing_file <- file.path(sample_info_dir, "disease-timing.tsv")

# path to sample metadata 
sample_metadata_file <- here::here("s3_files", "scpca-sample-metadata.tsv")

# color palette
diagnosis_group_palette <- here::here("palettes", "diagnosis-group-palette.tsv")

# output files 
pdf_dir <- here::here("figures", "pdfs")
output_pdf_file <- file.path(pdf_dir, "Fig1_sample-summary.pdf")


# Prep sample metadata ------------------------------------------------------

# read in project whitelist and grouping metadata 
project_whitelist <- readLines(project_whitelist_file)
diagnosis_groupings_df <- readr::read_tsv(diagnosis_groupings_file) |>
  dplyr::select(submitted_diagnosis, diagnosis_group)
disease_timing_df <- readr::read_tsv(disease_timing_file)

# read in sample metadata and filter to only projects in whitelist
sample_metadata_df <- readr::read_tsv(sample_metadata_file) |> 
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


diagnosis_count <- diagnosis_plot_df |> 
  dplyr::count(diagnosis) 

# create dataframe with diagnosis, diagnosis count, and disease timing for plotting
plot_df <- diagnosis_plot_df |> 
  dplyr::left_join(diagnosis_count) |>
  dplyr::mutate(diagnosis = forcats::fct_reorder(diagnosis, n),
                # relevel so other is last
                diagnosis_group = forcats::fct_relevel(diagnosis_group, 
                                                       "Other solid tumors",
                                                       after = Inf))
# Plot -------------------------------------------------------------------------

# get list of colors to use 
diagnosis_group_palette_df <- readr::read_tsv(diagnosis_group_palette)


# get list of all groups 
diagnosis_colors <- diagnosis_group_palette_df$color |> 
  purrr::set_names(diagnosis_group_palette_df$diagnosis_group)

# create faceted plot, one panel for each diagnosis group
# diagnosis is on the y axis and number of samples on the x axis
# bars are patterned based on the disease timing group 
diagnosis_plot <- ggplot(plot_df, aes(y = diagnosis,  fill = diagnosis_group)) +
  ggpattern::geom_bar_pattern(aes(pattern = forcats::fct_rev(standardized_disease_timing)), 
                              color = "black",
                              pattern_color = "black",
                              pattern_fill = "black",
                              pattern_density = 0.2,
                              pattern_spacing = 0.02)+
  facet_wrap(facets = "diagnosis_group", scales = "free") +
  # add label for the number of samples in each diagnosis group 
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = -0.25, size = 4) +
  # manually set patterns 
  ggpattern::scale_pattern_manual(values = c(
    "Initial diagnosis" = 'none', 
    "Post-mortem" = 'stripe', 
    "Progressive" = 'crosshatch', 
    "Recurrence" = 'circle',
    "Unknown" = 'wave'), 
    drop = FALSE) +
  # set colors to use palette for diagnosis groups
  scale_fill_manual(values = diagnosis_colors) +
  # make legend black and white 
  guides(fill = "none", 
         pattern = guide_legend(override.aes = list(fill = "white"), 
                                title.position = "top", 
                                title.hjust = 0.5,
                                reverse = TRUE)) + 
  labs(pattern = "Disease timing",
       x = "Number of samples",
       y = "Diagnosis") + 
  theme_classic() + 
  theme(text = element_text(size = 14),
        legend.position = "top",
        legend.key.size = unit(1, 'cm'),
        plot.margin = margin(1, 1, 1, 1, 'cm')
        
  ) +
  # make sure bar labels don't get cut off
  coord_cartesian(clip = "off")

# save plot 
ggsave(output_pdf_file, plot = diagnosis_plot, width = 15, height = 10)
