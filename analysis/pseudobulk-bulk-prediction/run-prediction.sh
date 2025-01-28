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
result_dir="results"

mkdir -p $scpca_dir
mkdir -p $tpm_dir
mkdir -p $result_dir

# Step 0: Sync data files from S3
Rscript ${script_dir}/sync-data-files.R \
  --output_dir ${scpca_dir}

for project_dir in $scpca_dir/*; do
    project_id=$(basename $project_dir)

    tpm_file="${tpm_dir}/${project_id}-tpm.tsv"

    # Step 1: Calculate TPM for each project
    Rscript ${script_dir}/calculate-tpm.R \
        --input_dir "${scpca_dir}/${project_id}" \
        --output_file "${tpm_file}"
done
