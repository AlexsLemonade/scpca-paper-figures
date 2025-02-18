#!/bin/bash

# This script runs an analysis using ESTIMATE

set -euo pipefail

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd "$basedir"

# Define directories
data_dir="data"
script_dir="scripts"
scpca_dir="${data_dir}/scpca"

mkdir -p $scpca_dir

ensembl_symbol_map_file="${data_dir}/ensembl-symbol-map.tsv"

# Sync data files from S3
Rscript ${script_dir}/sync-data-files.R \
  --output_dir "${scpca_dir}" \
  --ensembl_symbol_map_file "${ensembl_symbol_map_file}"

for project_dir in $scpca_dir/*; do

    project_id=$(basename $project_dir)
    expression_file="${data_dir}/${project_id}_expression.rds"

    # Calculate all expression quantities 
    Rscript ${script_dir}/prepare-expression-data.R \
      --project_id "${project_id}" \
      --scpca_dir "${project_dir}" \
      --output_file "${expression_file}"
done
