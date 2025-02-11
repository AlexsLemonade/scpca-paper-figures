#!/bin/bash

# This script runs the analysis predicting bulk from pseudobulk

set -euo pipefail

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd "$basedir"

# Define directories
data_dir="data"
script_dir="scripts"
scpca_dir="${data_dir}/scpca_data"
tpm_dir="${data_dir}/tpm"
pseudobulk_dir="${data_dir}/pseudobulk"
result_dir="results"
model_notebook_dir="model-notebooks"

mkdir -p $scpca_dir
mkdir -p $tpm_dir
mkdir -p $pseudobulk_dir
mkdir -p $result_dir

map_file="${data_dir}/bulk-library-sample-ids.tsv"


# Sync data files from S3
Rscript ${script_dir}/sync-data-files.R \
  --output_dir "${scpca_dir}" \
  --map_file "${map_file}"

# Prepare bulk counts data for comparisons
Rscript ${script_dir}/prepare-bulk-counts.R \
  --input_dir "${scpca_dir}" \
  --map_file "${map_file}" \
  --output_counts_file "${data_dir}/normalized-bulk-counts.rds" \
  --output_frac_expressed_file "${data_dir}/fraction-expressed-bulk.tsv"

for project_dir in $scpca_dir/*; do
    project_id=$(basename $project_dir)

    pseudobulk_file="${pseudobulk_dir}/${project_id}_pseudobulk.tsv"
    fraction_expressed_file="${data_dir}/${project_id}_fraction-expressed-single-cell.tsv"

    ###### TPMs are not currently used in the analysis ######
    # Calculate bulk TPM for each project
    #tpm_file="${tpm_dir}/${project_id}_tpm.tsv"
    #Rscript ${script_dir}/calculate-tpm.R \
    #  --input_dir "${scpca_dir}/${project_id}" \
    #  --output_pseudobulk_file "${tpm_file}"

    # Calculate pseudobulk matrices for each project
    Rscript ${script_dir}/calculate-pseudobulk.R \
      --input_dir "${scpca_dir}/${project_id}" \
      --output_pseudobulk_file "${pseudobulk_file}" \
      --output_frac_expressed_file "${fraction_expressed_file}"
done

# Build and export models to results/models with and without filtering unexpressed genes
for filter_genes in 0 1; do

  if [[ ${filter_genes} -eq 0 ]]; then
    filter_str="unfiltered"
  else
    filter_str="filtered"
  fi

  Rscript -e "rmarkdown::render('${model_notebook_dir}/build-assess-models.Rmd',
              params = list(filter_genes = ${filter_genes}),
              output_file = 'build-assess-models_${filter_str}.nb.html',
              output_dir = '${model_notebook_dir}')"
done