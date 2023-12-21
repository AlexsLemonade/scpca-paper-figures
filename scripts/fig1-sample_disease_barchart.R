# Generates a bar plot summarizing the diagnoses type in the ScPCA Portal 
# The resulting plot is faceted by the diagnoses group 
# Individual bars have different patterns based on the disease timing 

# Set up -----------------------------------------------------------------------

# load any libaries 
library(ggplot2)

# set up paths 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
plots_dir <- file.path(root_dir, "figures", "pngs")
output_plot_file <- file.path(plots_dir, "Fig1-sample-summary.png")

sample_metadata_s3 <- 's3://ccdl-scpca-data/sample_info/scpca-sample-metadata.tsv'
local_s3_dir <- file.path(root_dir, "s3_files")
fs::dir_create(local_s3_dir)

sample_info_dir <- file.path(root_dir, "sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
diagnosis_groupings_file <- file.path(sample_info_dir, "diagnosis-groupings.tsv")
disease_timing_file <- file.path(sample_info_dir, "disease-timing.tsv")

# read in metadata 
project_whitelist <- readLines(project_whitelist_file)
diagnosis_groupings_df <- readr::read_tsv(diagnosis_groupings_file) |>
  dplyr::select(submitted_diagnosis, diagnosis_group)
disease_timing_df <- readr::read_tsv(disease_timing_file)

# Copy and filter sample metadata 
local_sample_metadata <- file.path(local_s3_dir, "scpca-sample-metadata.tsv")

## TODO: revisit why this isn't actually working? 
sync_call <- paste('op run -- aws s3 cp', sample_metadata_s3, local_sample_metadata, sep = " ")
system(sync_call)

sample_metadata_df <- readr::read_tsv(local_sample_metadata) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist)

# Join sample metadata with diagnosis grouping and disease timing 
diagnosis_plot_df <- sample_metadata_df |> 
  dplyr::select(scpca_sample_id, diagnosis, disease_timing) |> 
  dplyr::left_join(diagnosis_groupings_df, by = c("diagnosis" = "submitted_diagnosis")) |> 
  dplyr::left_join(disease_timing_df, by = c("disease_timing" = "submitted_disease_timing"))

# get list of all groups 
diagnosis_groups <- diagnosis_plot_df |> 
  # right now let's remove this group 
  dplyr::filter(diagnosis_group != "Non-cancerous") |> 
  dplyr::pull(diagnosis_group) |>
  unique() |>
  forcats::fct_relevel("Brain and CNS",
                       "Leukemia", 
                       "Sarcoma",
                       "Other solid tumors")

# colors 
#names(diagnosis_groups) <- c("#201158", "#005E5E", "#578B21", "#E89E6B")
names(diagnosis_groups) <- c("#DF536B", "#61D04F", "#2297E6", "#28E2E5")


# create plot panels
barplot_panel <- function(diagnosis_plot_df, diagnosis_group, group_color){
  
  group_to_plot <- diagnosis_plot_df |> 
    dplyr::filter(diagnosis_group %in% {{diagnosis_group}}) 
  
  diagnosis_count <- group_to_plot |> 
    dplyr::count(diagnosis) 
  
  plot_df <- group_to_plot |> 
    dplyr::left_join(diagnosis_count) |>
    dplyr::mutate(diagnosis = forcats::fct_reorder(diagnosis, -n))
    
  
  diagnosis_plot <- ggplot(plot_df, aes(y = diagnosis)) + 
    ggpattern::geom_bar_pattern(aes(pattern = standardized_disease_timing),
                                pattern_density = 0.04,
                                fill = group_color,
                                color = "black",
                                pattern_color = "black",
                                pattern_angle = 0,
                                pattern_spacing = 0.01) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90),
          text = element_text(size = 10),
          legend.key.size = unit(1, 'cm'),
          legend.key = element_blank(),
          legend.position = 'none') +
    labs(
      x = "",
      y = "",
      title = diagnosis_group,
    ) +
    ggpattern::scale_pattern_manual(values = c('circle', 'stripe', 'wave', 'crosshatch', 'plain'))

  return(diagnosis_plot)
}

all_panels <- diagnosis_groups |> 
  purrr::imap(\(group, group_color){
    barplot_panel(diagnosis_plot_df = diagnosis_plot_df,
                  diagnosis_group = group, 
                  group_color = group_color)
  })

combined_plot <- patchwork::wrap_plots(all_panels) +
  patchwork::plot_annotation(tag_levels = 'A')
  

ggsave(output_plot_file, width = 15, height = 10)
 

# combine with patchwork 

facet_plot <- ggplot(diagnosis_plot_df, aes(y = diagnosis, fill = diagnosis_group)) + 
  ggpattern::geom_bar_pattern(aes(pattern = standardized_disease_timing),
                              pattern_density = 0.04,
                              color = "black",
                              pattern_color = "black",
                              pattern_angle = 0,
                              pattern_spacing = 0.01) +
  theme_classic() + 
  facet_wrap(vars(diagnosis_group)) +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        legend.key.size = unit(1, 'cm')) +
  labs(
    x = "",
    y = ""
  ) +
  ggpattern::scale_pattern_manual(values = c('circle', 'stripe', 'wave', 'crosshatch', 'circle', 'weave'))

ggsave("faceted_plot.png")

