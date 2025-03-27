#!/usr/bin/env Rscript

# This script is used to create a table with the brain subdiagnoses (HGG vs. LGG)
# for all samples included in Fig5B and 5C 
# The output is a table with diagnosis and subdiagnosis_group saved to `sample-info/brain-classifications-no-multiplexed.tsv`

# SCPCP000001 samples will be labeled as HGG
# SCPCP000002 samples will be labeled as LGG
# multiplexed samples are not included
# all other samples are assigned based on the category in https://pmc.ncbi.nlm.nih.gov/articles/PMC8328013/#T1

# define input metadata
sample_metadata_file <- here::here("s3_files", "scpca-sample-metadata.tsv")
library_metadata_file <- here::here("s3_files", "scpca-library-metadata.tsv")

# define output file 
classification_output_file <- here::here("sample-info", "brain-classifications-no-multiplexed.tsv")

# brain projects 
brain_project_ids <- c("SCPCP000001", "SCPCP000002", "SCPCP000010", "SCPCP000021", "SCPCP000009")
# pull out those that are non-multiplex single cell/nuc
non_multiplex_samples <- readr::read_tsv(library_metadata_file) |> 
  dplyr::filter(!stringr::str_detect(scpca_sample_id, ";"),
                seq_unit %in% c("cell", "nucleus")) |> 
  dplyr::pull(scpca_sample_id)

# read in sample metadata and select samples
sample_df <- readr::read_tsv(sample_metadata_file) |> 
  dplyr::filter(scpca_project_id %in% brain_project_ids & scpca_sample_id %in% non_multiplex_samples)

# define hgg vs lgg types
hgg_types <- c(
  "Anaplastic glioma", # part of SCPCP000001
  "Glioblastoma", # part of SCPCP000001
  "High-grade glioma", # part of SCPCP000001
  "Diffuse midline glioma",# part of SCPCP000001
  "Pleomorphic xanthoastrocytoma", # part of SCPCP000001
  "Anaplastic astrocytoma", # part of SCPCP000001
  "Infant-type hemispheric glioma"
)

lgg_types <- c(
  "Pilocytic astrocytoma", # part of SCPCP000002
  "Ganglioglioma", # part of SCPCP000002
  "Low-grade glioma", # part of SCPCP000002
  "Ganglioglioma/ATRT", # part of SCPCP000002
  "Ependymoma", # part of SCPCP000002
  "Dysplastic gangliocytoma", 
  "Dysembryoplastic neuroepithelial tumor", 
  "Glial-neuronal tumor"
)

# assign classifications
classification_df <- sample_df |> 
  # get rid of the normal sample
  dplyr::filter(diagnosis != "Non-cancerous") |> 
  # assign to HGG vs LGG 
  dplyr::mutate(
    subdiagnosis_group = dplyr::case_when(
      diagnosis %in% hgg_types ~ "High-grade glioma",
      diagnosis %in% lgg_types ~ "Low-grade glioma",
      .default = "other"
    )
  ) |> 
  dplyr::select(diagnosis, subdiagnosis_group) |> 
  unique() 

readr::write_tsv(classification_df, classification_output_file)

