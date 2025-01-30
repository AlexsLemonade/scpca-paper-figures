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

mkdir -p $scpca_dir
mkdir -p $tpm_dir
mkdir -p $pseudobulk_dir
mkdir -p $result_dir

# Step 0: Sync data files from S3
Rscript ${script_dir}/sync-data-files.R \
  --output_dir ${scpca_dir}

for project_dir in $scpca_dir/*; do
    project_id=$(basename $project_dir)

    tpm_file="${tpm_dir}/${project_id}-tpm.tsv"
    pseudobulk_file="${pseudobulk_dir}/${project_id}-pseudobulk.tsv"

    # Step 1: Calculate bulk TPM for each project
    if [ ! -f ${tpm_file} ]; then
      Rscript ${script_dir}/calculate-tpm.R \
        --input_dir "${scpca_dir}/${project_id}" \
        --output_file "${tpm_file}"
    fi

    # Step 2: Calculate pseudobulk matrices for each project
    if [ ! -f ${pseudobulk_file} ]; then
      Rscript ${script_dir}/calculate-pseudobulk.R \
        --input_dir "${scpca_dir}/${project_id}" \
        --output_file "${pseudobulk_file}"
    fi
done
