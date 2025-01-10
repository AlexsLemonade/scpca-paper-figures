# Generates a bar plot summarizing the diagnoses type in the ScPCA Portal 
# The resulting plot is faceted by the diagnoses group 
# Individual bars have different patterns based on the disease timing 

renv::load()

# load any libaries 
library(ggplot2)

# Set up -----------------------------------------------------------------------


# source in helper functions for plotting
function_file <- here::here("scripts", "utils", "sample-summary-helper-functions.R")
source(function_file)

# all metadata files 
sample_info_dir <- here::here("sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
diagnosis_groupings_file <- file.path(sample_info_dir, "diagnosis-groupings.tsv")
disease_timing_file <- file.path(sample_info_dir, "disease-timing.tsv")

# path to metadata files
sample_metadata_file <- here::here("s3_files", "scpca-sample-metadata.tsv")
project_metadata_file <- here::here("s3_files", "scpca-project-metadata.tsv")

# color palette
disease_timing_palette <- here::here("palettes", "disease-timing-palette.tsv")

# output files 
pdf_dir <- here::here("figures", "pdfs")
output_pdf_file <- file.path(pdf_dir, "Fig1_sample-summary.pdf")

tables_dir <- here::here("manuscript-numbers")
diagnosis_count_table <- file.path(tables_dir, "diagnosis-group-counts.tsv")
disease_timing_count_table <- file.path(tables_dir, "disease-timing-counts.tsv")

# Prep sample metadata ------------------------------------------------------

# read in project whitelist and grouping metadata 
project_whitelist <- readLines(project_whitelist_file)

# get sample whitelist 
sample_whitelist <- get_sample_whitelist(project_metadata_file, project_whitelist)

# read in diagnosis groupings and disease timing
diagnosis_groupings_df <- readr::read_tsv(diagnosis_groupings_file) |>
  dplyr::select(submitted_diagnosis, diagnosis_group)
disease_timing_df <- readr::read_tsv(disease_timing_file)

# read in sample metadata and filter to only projects in whitelist
sample_metadata_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist, 
                scpca_sample_id %in% sample_whitelist,
                # remove any cell lines from this plot
                !is_cell_line)

# Join sample metadata with diagnosis grouping and disease timing 
diagnosis_plot_df <- sample_metadata_df |> 
  dplyr::select(scpca_sample_id, diagnosis, disease_timing) |> 
  dplyr::left_join(diagnosis_groupings_df, by = c("diagnosis" = "submitted_diagnosis")) |> 
  dplyr::left_join(disease_timing_df, by = c("disease_timing" = "submitted_disease_timing")) |> 
  # right now let's remove this group
  dplyr::filter(diagnosis_group != "Non-cancerous") |> 
  # make it a factor so we can use ggpattern 
  dplyr::mutate(standardized_disease_timing = as.factor(standardized_disease_timing) |> 
                  forcats::fct_rev() |> 
                  forcats::fct_relevel("During or after treatment", after = 1))


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

# tables counting total for each group and disease timing
diagnosis_group_totals <- plot_df |> 
  dplyr::count(diagnosis_group)
readr::write_tsv(diagnosis_group_totals, diagnosis_count_table)

disease_timing_totals <- plot_df |> 
  dplyr::count(standardized_disease_timing)
readr::write_tsv(disease_timing_totals, disease_timing_count_table)
# Plot -------------------------------------------------------------------------

# get list of colors to use 
disease_timing_palette_df <- readr::read_tsv(disease_timing_palette)


# get list of all groups 
disease_timing_colors <- disease_timing_palette_df$color |> 
  purrr::set_names(disease_timing_palette_df$disease_timing)

# create faceted plot, one panel for each diagnosis group
# diagnosis is on the y axis and number of samples on the x axis
# bars are colored based on the disease timing group 
diagnosis_plot <- ggplot(plot_df, aes(y = diagnosis,  fill = disease_timing)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(facets = "diagnosis_group", scales = "free") +
  # add label for the number of samples in each diagnosis group 
  geom_text(aes(label = after_stat(count)), stat = "count", hjust = -0.25, size = 4) +
  # set colors to use palette for diagnosis groups
  scale_fill_manual(values = disease_timing_colors) +
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
  theme(text = element_text(size = 12),
        legend.position = "top",
        legend.key.size = unit(1, 'cm'),
        plot.margin = margin(1, 1, 1, 1, 'cm')
        
  ) +
  # make sure bar labels don't get cut off
  coord_cartesian(clip = "off")

# save plot 
ggsave(output_pdf_file, plot = diagnosis_plot, width = 15, height = 10)
