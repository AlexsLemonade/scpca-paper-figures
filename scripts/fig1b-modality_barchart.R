# Generates a bar plot summarizing the modalities represented in all libraries
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

# TODO: create a color palette for single-cell/single-nuclei (?)

# output files 
plots_dir <- file.path(root_dir, "figures", "pngs")
output_plot_file <- file.path(plots_dir, "Fig1B-modality-summary.png")

# Prep sample metadata ------------------------------------------------------

# read in project whitelist
project_whitelist <- readLines(project_whitelist_file)

# read in library metadata and filter to only projects in whitelist
library_metadata_df <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist) |>
  dplyr::select(scpca_project_id, scpca_sample_id, scpca_library_id, seq_unit, technology)

# create a modality column that labels everything as single suspension, bulk, spatial, or CITE
modality_df <- library_metadata_df |> 
  dplyr::mutate(
    modality = dplyr::case_when(
      seq_unit %in% c("cell", "nucleus")  & !(stringr::str_detect(technology,"CITE")) ~ "Single suspension",
      seq_unit == "bulk" ~ "bulk",
      seq_unit == "spot" ~ "spatial",
      stringr::str_detect(technology,"CITE") ~ "CITE"
    )
  )

# get a list of any bulk samples that don't have matching single suspension 
all_bulk <- modality_df |> 
  dplyr::filter(modality == "bulk") |> 
  dplyr::pull(scpca_sample_id)

all_spatial <- modality_df |> 
  dplyr::filter(modality == "spatial") |>
  dplyr::pull(scpca_sample_id)

# get list of all single cell and all single nuclei
all_single_cell <- modality_df |> 
  dplyr::filter(seq_unit == "cell" & modality == "Single suspension") |> 
  dplyr::pull(scpca_sample_id)
all_single_nuc <- modality_df |> 
  dplyr::filter(seq_unit == "nucleus" & modality == "Single suspension") |> 
  dplyr::pull(scpca_sample_id)

# get list of all bulk and spatial samples that don't have a corresponding cell/nucleus sample
all_suspension <- c(all_single_cell, all_single_nuc)
bulk_only <- setdiff(all_bulk, all_suspension)
spatial_only <- setdiff(all_spatial, c(all_single_cell, all_single_nuc))


# remove bulk and spatial only from metadata
filtered_modality_df <- modality_df |> 
  dplyr::filter(!(scpca_sample_id %in% c(bulk_only, spatial_only))) |> 
  # make sure all bulk and spatial get set with cell or nucleus
  dplyr::mutate(
    seq_unit = dplyr::case_when(
      scpca_sample_id %in% all_single_cell ~ "cell",
      scpca_sample_id %in% all_single_nuc ~ "nucleus"
    )
  ) |>
  dplyr::group_by(modality) |> 
  dplyr::add_tally(name = "total_per_modality")

# TODO: How do we turn in the number of libraries to a percentage 
# based on the percentage of total single-cell/single-nuc libraries

# make plot 
ggplot(filtered_modality_df, aes(x = forcats::fct_reorder(modality, -total_per_modality), fill = seq_unit)) +
  geom_bar(stat = "count") +
  theme_classic() + 
  labs(
    x = "Type of library",
    y = "Number of libraries",
    fill = ""
  )
  