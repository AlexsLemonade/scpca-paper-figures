# This script creates a table showing the compilation of all references 
# built for CellAssign from PanglaoDB
# specifically one row per project, with list of diagnoses, ref name, and organ list

renv::load()

# Set up -----------------------------------------------------------------------

# sample metadata file
# need to summarize diagnoses for each project
sample_metadata_file <- here::here("s3_files", "scpca-sample-metadata.tsv")

# panglao reference metadata
# need this to get list of organs and ref name
panglao_metadata_file <- here::here("sample-info", "celltype-reference-metadata.tsv")

# project celltype metadata 
# need this to get ref name for each project
project_metadata_file <- here::here("s3_files", "scpca-project-celltype-metadata.tsv")

# output table 
table_dir <- here::here("tables")
fs::dir_create(table_dir)
output_table_file <- file.path(table_dir, "TableS2_cellassign-ref-summary.tsv")

# Prepare table ----------------------------------------------------------------

# read in all metadata 
sample_df <- readr::read_tsv(sample_metadata_file)
project_df <- readr::read_tsv(project_metadata_file)
# only keep ref name and organs for joining later 
panglao_df <- readr::read_tsv(panglao_metadata_file) |> 
  dplyr::select(celltype_ref_name, organs)

# filter project df to include only projects with a cellassign ref
project_df <- project_df |> 
  dplyr::filter(!is.na(cellassign_ref_name)) |> 
  # only keep project ID and ref name for joining 
  dplyr::select(scpca_project_id, cellassign_ref_name)

grouped_sample_df <- sample_df |> 
  dplyr::select(scpca_project_id, diagnosis) |> 
  # remove any projects that don't have cell assign refs
  dplyr::filter(scpca_project_id %in% project_df$scpca_project_id) |> 
  dplyr::group_by(scpca_project_id) |> 
  # create a list of diagnoses for each project
  dplyr::summarize(Diagnoses = paste(unique(diagnosis), collapse = ", ")) |> 
  # add ref name
  dplyr::left_join(project_df) |> 
  # add organs
  dplyr::left_join(panglao_df, by = c("cellassign_ref_name" = "celltype_ref_name")) |> 
  dplyr::rename(
    "ScPCA Project ID" = "scpca_project_id",
    "PanglaoDB reference name" = "cellassign_ref_name",
    "Organs included in reference" = "organs"
  )

# export table 
readr::write_tsv(grouped_sample_df, output_table_file)

