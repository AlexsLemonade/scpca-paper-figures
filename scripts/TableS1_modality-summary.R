# This script creates a table summarizing the total number of libraries and types 
# of libraries for each project 

renv::load()

# load any libaries 
library(ggplot2)

# Set up -----------------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# path to metadata files
s3_file_dir <- file.path(root_dir, "s3_files")
library_metadata_file <- file.path(s3_file_dir, "scpca-library-metadata.tsv")
sample_metadata_file <- file.path(s3_file_dir, "scpca-sample-metadata.tsv")
project_metadata_file <- file.path(s3_file_dir, "scpca-project-metadata.tsv")

# project whitelist and diagnosis groupings
sample_info_dir <- file.path(root_dir, "sample-info")
project_whitelist_file <- file.path(sample_info_dir, "project-whitelist.txt")
diagnosis_groupings_file <- file.path(sample_info_dir, "diagnosis-groupings.tsv")

# output files 
table_dir <- file.path(root_dir, "tables")
fs::dir_create(table_dir)
output_table_file <- file.path(table_dir, "TableS1-modality-overview.tsv")

# read in project whitelist
project_whitelist <- readLines(project_whitelist_file)

# read in project metadata files and create sample whitelist 
project_metadata_files <- readr::read_tsv(project_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist) |> 
  dplyr::pull(metadata_file)
# get full file paths to each project metadata file
project_metadata_files <- file.path(s3_file_dir, "project-metadata", project_metadata_files)

# grab samples that are on the portal and create a whitelist 
sample_whitelist <- project_metadata_files |> 
  purrr::map(\(file){
    sample_list <- readr::read_tsv(file) |> 
      dplyr::filter(on_portal) |> 
      dplyr::pull(scpca_sample_id)
    
    return(sample_list)
  }) |> 
  unlist()

# read in groupings
diagnosis_groupings_df <- readr::read_tsv(diagnosis_groupings_file)

# read in library metadata and filter to included projects
library_metadata_df <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist)

# read in sample metadata and filter to included projects and samples
sample_metadata_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% project_whitelist,
                scpca_sample_id %in% sample_whitelist) |> 
  dplyr::select(scpca_project_id, scpca_sample_id, diagnosis) |> 
  # add diagnosis groupings
  dplyr::left_join(diagnosis_groupings_df, by = c("diagnosis" = "submitted_diagnosis"))

# Create table -----------------------------------------------------------------

demuxed_metadata_df <- library_metadata_df |> 
  # make sure we have one row per library/ sample combination
  # this specifically separates sample IDs for multiplex libraries for easier joining to sample metadata
  dplyr::mutate(scpca_sample_id = stringr::str_split(scpca_sample_id, ";")) |> 
  tidyr::unnest(scpca_sample_id) |> 
  # now that sample ids are separated, filter to whitelist 
  dplyr::filter(scpca_sample_id %in% sample_whitelist)
  
# create a total count of the number of libraries per project per diagnosis group
total_library_count <- demuxed_metadata_df |>
  dplyr::left_join(sample_metadata_df) |> 
  dplyr::select(scpca_project_id, scpca_library_id, diagnosis_group) |> 
  # remove any duplicates (e.g. those with CITE or with cell hash)
  dplyr::distinct() |>
  dplyr::count(scpca_project_id, diagnosis_group, name = "Total number of libraries")

total_sample_count <- demuxed_metadata_df |>
  dplyr::left_join(sample_metadata_df) |> 
  dplyr::select(scpca_project_id, scpca_sample_id, diagnosis_group) |> 
  # remove any duplicates (e.g. those with CITE or with cell hash)
  dplyr::distinct() |> 
  dplyr::count(scpca_project_id, diagnosis_group, name = "Total number of samples")

modality_counts_df <- demuxed_metadata_df |> 
  dplyr::mutate(
    # make a new column that summarizes all modality info
    # options are cell, nucleus, bulk, spatial, cite, and multiplexed
    # the assigned values will eventually be the column names used in the output
    modality = dplyr::case_when(
      seq_unit == "bulk" ~ "Bulk RNA",
      seq_unit == "spot" ~ "Spatial transcriptomics",
      stringr::str_detect(technology,"CITE") ~ "With CITE-seq",
      stringr::str_detect(technology,"cellhash") ~ "With cell hashing",
      seq_unit == "cell" ~ "Single-cell",
      seq_unit == "nucleus" ~ "Single-nucleus"
    )
  ) |> 
  # join with sample metadata to add in diagnosis groupings
  dplyr::left_join(sample_metadata_df) |>
  dplyr::select(
    scpca_project_id, scpca_library_id, modality, diagnosis, diagnosis_group
  ) 

# determine full diagnoses for each project
full_diagnoses <- modality_counts_df |>
  # for each project, group by diagnosis group, and summarize the diagnoses included
  dplyr::group_by(scpca_project_id, diagnosis_group) |> 
  dplyr::summarize(diagnosis = paste(unique(diagnosis), collapse = ";"))

# get total numbers of libraries for each 10X kit for each project 
kit_counts <- demuxed_metadata_df |> 
  dplyr::left_join(sample_metadata_df) |>
  dplyr::select(scpca_project_id, scpca_library_id, diagnosis_group, technology) |> 
  dplyr::distinct() |> 
  dplyr::mutate(    # create a separate column with the kit used for each sample
    kit_type = dplyr::case_when(
      # account for both 10Xv2 3' and 5'
      technology %in% c("10Xv2", "10Xv2_5prime") ~ "10Xv2",
      technology == "10Xv3" ~ "10Xv3",
      technology == "10Xv3.1" ~ "10Xv3.1"
    )) |> 
  # remove any NA, these are cell hash, bulk, ST, and CITE
  tidyr::drop_na(kit_type) |> 
  # count each combination of project, diagnosis group and kit
  dplyr::count(scpca_project_id, diagnosis_group, kit_type) |> 
  tidyr::pivot_wider(names_from = kit_type, values_from = n, values_fill = 0)

summarized_counts_df <- modality_counts_df |>
  dplyr::select(scpca_project_id, scpca_library_id, diagnosis_group, modality) |>
  dplyr::distinct() |>
  # count each combination of diagnosis group and modality
  dplyr::count(scpca_project_id, diagnosis_group, modality) |>
  # create a column for each modality that contains the total number of libraries
  tidyr::pivot_wider(names_from = modality, values_from = n, values_fill = 0) |>
  dplyr::left_join(total_library_count) |>
  dplyr::left_join(total_sample_count) |> 
  dplyr::left_join(full_diagnoses) |>
  dplyr::left_join(kit_counts) |> 
  # set desired order and do some renaming
  dplyr::relocate(
    "Diagnosis group" = diagnosis_group, 
    "Diagnoses" = diagnosis, 
    "Total number of samples (S)" = "Total number of samples",
    "Total number of libraries (L)" = "Total number of libraries", 
    "Single-cell (L)" = "Single-cell", 
    "Single-nucleus (L)" = "Single-nucleus", 
    "10Xv2 (L)" = "10Xv2",
    "10Xv3 (L)" = "10Xv3",
    "10Xv3.1 (L)" = "10Xv3.1",
    "Bulk RNA (L)" = "Bulk RNA",
    "Spatial transcriptomics (L)" = "Spatial transcriptomics", 
    "With CITE-seq (L)" = "With CITE-seq",
    "With cell hashing (L)" = "With cell hashing",
    .after = 1
    ) 


# export table 
readr::write_tsv(summarized_counts_df, output_table_file)

