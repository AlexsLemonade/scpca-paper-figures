# Generates a bar plot summarizing the modalities represented in all samples
# found on the portal 

# load any libaries 
library(ggplot2)

# Set up -----------------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# path to metadata files
library_metadata_file <- file.path(root_dir, "s3_files", "scpca-library-metadata.tsv")
sample_metadata_file <- file.path(root_dir, "s3_files", "scpca-sample-metadata.tsv")

# project whitelist file 
project_whitelist_file <- file.path(root_dir, "sample-info", "project-whitelist.txt")

# color palette for single-cell/single-nuclei
suspension_palette_file <- file.path(root_dir, "palettes", "suspension-palette.tsv")

# output files 
plots_dir <- file.path(root_dir, "figures", "pngs")
output_plot_file <- file.path(plots_dir, "Fig1B-modality-summary.png")

# Prep metadata ------------------------------------------------------

# read in project whitelist
project_whitelist <- readLines(project_whitelist_file)

# read in color palette
suspension_palette <- readr::read_tsv(suspension_palette_file)

# read in library metadata and filter to only projects in whitelist
library_metadata_df <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist) |>
  dplyr::select(scpca_project_id, scpca_sample_id, scpca_library_id, seq_unit, technology) |> 
  # create a modality column that labels everything as single suspension, bulk, spatial, or CITE
  dplyr::mutate(
    modality = dplyr::case_when(
      # all libraries that are either cell or nucleus but not CITE or hash
      seq_unit %in% c("cell", "nucleus")  & 
        (!(stringr::str_detect(technology,"CITE")) & 
            !(stringr::str_detect(technology, "cellhash"))) ~ "Single suspension",
      seq_unit == "bulk" ~ "Bulk",
      seq_unit == "spot" ~ "Spatial transcriptomics",
      stringr::str_detect(technology,"CITE") ~ "With CITE-seq",
      stringr::str_detect(technology, "cellhash") ~ "With cell hashing"
    )
  )

# get a list of all bulk and spatial samples
# these will be checked against all single-cell/ nuc samples to look for matching samples
all_bulk <- library_metadata_df |> 
  dplyr::filter(modality == "Bulk") |> 
  dplyr::pull(scpca_sample_id)

all_spatial <- library_metadata_df |> 
  dplyr::filter(modality == "Spatial transcriptomics") |>
  dplyr::pull(scpca_sample_id)

# get list of all single cell and all single nuclei
all_single_cell <- library_metadata_df |> 
  dplyr::filter(seq_unit == "cell" & modality == "Single suspension") |> 
  dplyr::pull(scpca_sample_id)
all_single_nuc <- library_metadata_df |> 
  dplyr::filter(seq_unit == "nucleus" & modality == "Single suspension") |> 
  dplyr::pull(scpca_sample_id)

# get list of all bulk and spatial samples that don't have a corresponding cell/nucleus sample
all_suspension <- c(all_single_cell, all_single_nuc)
bulk_only <- setdiff(all_bulk, all_suspension)
spatial_only <- setdiff(all_spatial, all_suspension)


# remove bulk and spatial only from metadata
filtered_modality_df <- library_metadata_df |> 
  dplyr::filter(!(scpca_sample_id %in% c(bulk_only, spatial_only))) |>
  # make sure all bulk and spatial get designated with cell or nucleus
  # also rename to be more specific when creating legend
  dplyr::mutate(
    seq_unit = dplyr::case_when(
      scpca_sample_id %in% all_single_cell ~ "Single-cell",
      scpca_sample_id %in% all_single_nuc ~ "Single-nuclei"
    )
  ) |> 
  # first combine all modalities for each sample id into one list
  # make sure to keep seq unit
  dplyr::group_by(scpca_sample_id, seq_unit) |>
  dplyr::summarise(modality = paste(modality, collapse = ";")) |> 
  # split any lists of modalities so we have one row per sample per modality
  dplyr::mutate(modality = stringr::str_split(modality, ";")) |> 
  tidyr::unnest(modality) |>
  dplyr::group_by(modality) |> 
  # get the total for each modality to use for specifying order in plot
  dplyr::add_tally(name = "total_per_modality") |>
  # add a column to help pull out additional modalities into its own facet 
  dplyr::mutate(additional_mods = ifelse(modality == "Single suspension", "All samples", "Samples with additional modalities"))

# Plot -------------------------------------------------------------------------

# get the colors for each suspension type
suspension_colors <- suspension_palette$color |> 
  purrr::set_names(suspension_palette$suspension_type)

# bar chart showing the total number of samples for each modality
# the bars are colored by what suspension type that modality is composed of
ggplot(filtered_modality_df, aes(x = forcats::fct_reorder(modality, -total_per_modality), fill = seq_unit)) +
  geom_bar(stat = "count") +
  # split the graph to designate all samples vs. those with additional modalities
  facet_grid(cols = vars(additional_mods),
             scales = "free",
             space = "free") +
  theme_classic() + 
  labs(
    x = "",
    y = "Number of samples",
    fill = "Suspension type"
  ) + 
  scale_fill_manual(values= suspension_colors) +
  theme(legend.position = c(.85,.9),
        legend.direction = "horizontal",
        legend.box.background = element_rect(colour = "black"),
        text = element_text(size = 14)) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave(output_plot_file)
  