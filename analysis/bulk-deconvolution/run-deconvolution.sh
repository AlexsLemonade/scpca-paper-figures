#!/bin/bash

# This script runs the bulk deconvolution analysis.

set -euo pipefail

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd "$basedir"

# Define directories
data_dir="data"
script_dir="scripts"
salmon_quant_dir="${data_dir}/salmon-quant-files"
reference_dir=${data_dir}/reference # stores the id map file
tpm_dir="${data_dir}/tpm"
result_dir="results"

mkdir -p $salmon_quant_dir
mkdir -p $reference_dir
mkdir -p $result_dir

# Step 0: Prepare id map and sync quant.sf files if they do not exist
Rscript ${script_dir}/prepare-data-files.R

for project_dir in $salmon_quant_dir/*; do
    project_id=$(basename $project_dir)

    tpm_file="${tpm_dir}/${project_id}-tpm.rds"
    quantiseq_file="${result_dir}/${project_id}-quantiseq.tsv"
    epic_file="${result_dir}/${project_id}-epic.tsv"

    # Step 1: Calculate TPM for each project
    Rscript ${script_dir}/calculate-tpm.R \
        --input_dir "${salmon_quant_dir}/${project_id}" \
        --output_file "${tpm_file}"

    # Step 2: Run quanTIseq on each project
    Rscript ${script_dir}/run-quantiseq.R \
        --input_file ${tpm_file} \
        --output_file ${quantiseq_file}

    # Step 3: Run EPIC on each project
    Rscript ${script_dir}/run-epic.R \
        --input_file ${tpm_file} \
        --output_file ${epic_file}
done
