# Generates a bar plot summarizing the modalities represented in all samples
# found on the portal 

renv::load()

# load any libaries 
library(ggplot2)

# Set up -----------------------------------------------------------------------

# path to metadata files
library_metadata_file <- here::here( "s3_files", "scpca-library-metadata.tsv")

# whitelist files 
project_whitelist_file <- here::here("sample-info", "project-whitelist.txt")
sample_whitelist_file <- here::here("sample-info", "sample-whitelist.txt")

# color palette for single-cell/single-nuclei
suspension_palette_file <- here::here("palettes", "suspension-palette.tsv")

# output files 
pdf_dir <- here::here("figures", "pdfs")
output_pdf_file <- file.path(pdf_dir, "Fig1B_modality-summary.pdf")

# set order of modalities for final plot 
modality_order <- c("Single suspension",
                    "Bulk",
                    "Spatial transcriptomics",
                    "With CITE-seq",
                    "With cell hashing"
                    )

# Prep metadata ------------------------------------------------------

# read in color palette
suspension_palette <- readr::read_tsv(suspension_palette_file)

# read in project and sample whitelists
project_whitelist <- readLines(project_whitelist_file)
sample_whitelist <- readLines(sample_whitelist_file)

# read in library metadata
library_metadata_df <- readr::read_tsv(library_metadata_file)


library_metadata_df <- library_metadata_df |> 
  dplyr::filter(scpca_project_id %in% project_whitelist) |>
  # separate out multiplexed samples
  dplyr::mutate(scpca_sample_id = stringr::str_split(scpca_sample_id, ";")) |> 
  tidyr::unnest(scpca_sample_id) |> 
  # now that sample ids are separated, filter to whitelist 
  dplyr::filter(scpca_sample_id %in% sample_whitelist) |> 
  # get one row per sample, seq unit, and tech combo 
  dplyr::select(scpca_project_id, scpca_sample_id, seq_unit, technology) |> 
  # create a modality column that labels everything as single suspension, bulk, spatial, or CITE
  dplyr::mutate(
    modality = dplyr::case_when(
      seq_unit == "bulk" ~ "Bulk",
      seq_unit == "spot" ~ "Spatial transcriptomics",
      stringr::str_detect(technology,"CITE") ~ "With CITE-seq",
      stringr::str_detect(technology, "cellhash") ~ "With cell hashing",
      seq_unit %in% c("cell", "nucleus") ~ "Single suspension"
    )
  ) |> 
  unique()


# get a list of all bulk and spatial samples
# these will be checked against all single-cell/ nuc samples to look for matching samples
all_bulk <- library_metadata_df |> 
  dplyr::filter(modality == "Bulk") |> 
  dplyr::pull(scpca_sample_id) |> 
  unique()

all_spatial <- library_metadata_df |> 
  dplyr::filter(modality == "Spatial transcriptomics") |>
  dplyr::pull(scpca_sample_id) |> 
  unique()

# get list of all single cell and all single nuclei
all_single_cell <- library_metadata_df |> 
  dplyr::filter(seq_unit == "cell" & modality == "Single suspension") |> 
  dplyr::pull(scpca_sample_id) |>
  unique()
all_single_nuc <- library_metadata_df |> 
  dplyr::filter(seq_unit == "nucleus" & modality == "Single suspension") |> 
  dplyr::pull(scpca_sample_id) |>
  unique()

# get list of all bulk and spatial samples that don't have a corresponding cell/nucleus sample
all_suspension <- c(all_single_cell, all_single_nuc) |> 
  unique() # account for two libraries that have both cell and nucleus
bulk_only <- setdiff(all_bulk, all_suspension)
spatial_only <- setdiff(all_spatial, all_suspension)


# remove bulk and spatial only from metadata
filtered_modality_df <- library_metadata_df |> 
  dplyr::filter(!(scpca_sample_id %in% c(bulk_only, spatial_only))) |>
  # two samples were sequenced using different 10x kits, so we don't want to double count them 
  dplyr::select(-technology) |> 
  # make sure all bulk and spatial get designated with cell or nucleus
  # also rename to be more specific when creating legend
  dplyr::mutate(
    seq_unit = dplyr::case_when(
      # two samples are sequenced using both cell and nuclei so count them separately 
      scpca_sample_id %in% c("SCPCS000250", "SCPCS000251") & seq_unit == "cell" ~ "Single-cell",
      scpca_sample_id %in% c("SCPCS000250", "SCPCS000251") & seq_unit == "nucleus" ~ "Single-nuclei",
      scpca_sample_id %in% all_single_cell ~ "Single-cell",
      scpca_sample_id %in% all_single_nuc ~ "Single-nuclei"
    )
  ) |> 
  unique() |> 
  # counts for seq unit and modality combos 
  dplyr::add_count(seq_unit, modality, name = "total_per_modality") |> 
  # add a column to help pull out additional modalities into its own facet 
  dplyr::mutate(
    additional_mods = ifelse(
      modality == "Single suspension", 
      "All samples", 
      "Samples with additional modalities"
    ),
    # set custom order
    # with cite and with cell hash should be together 
    modality = forcats::fct_relevel(modality, modality_order))


# Plot -------------------------------------------------------------------------

# get the colors for each suspension type
suspension_colors <- suspension_palette$color |> 
  purrr::set_names(suspension_palette$suspension_type)

# bar chart showing the total number of samples for each modality
# the bars are colored by what suspension type that modality is composed of
ggplot(filtered_modality_df, aes(x = modality, fill = seq_unit)) +
  geom_bar(stat = "count", color = "black") +
  # split the graph to designate all samples vs. those with additional modalities
  facet_grid(cols = vars(additional_mods),
             scales = "free",
             space = "free") +
  theme_classic() + 
  # add label to middle of each bar 
  geom_text(aes(label = after_stat(count)), 
            stat = "count", 
            position = position_stack(vjust = 0.5)) +
  labs(
    x = "",
    y = "Number of samples",
    fill = "Suspension type"
  ) + 
  scale_fill_manual(values= suspension_colors) +
  theme(legend.position.inside = c(.8,.88),
        legend.direction = "horizontal",
        legend.box.background = element_rect(color = "black"),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(position = "inside", title.position = "top", title.hjust = 0.5))

ggsave(output_pdf_file, width = 7, height = 7)
