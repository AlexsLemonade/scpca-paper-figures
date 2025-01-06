#!/bin/bash

# This script runs this analysis and I will add more comments later.
# If you are using 1P, you may need op run -- .

set -euo pipefail

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd "$basedir"

# Define directories
data_dir="data"
script_dir="scripts"
salmon_quant_dir="${data_dir}/salmon-quant-files"
tpm_dir="${data_dir}/tpm"

# Step 0: Download quant.sf files if they do not exist
if [ ! -d $salmon_quant_dir ]; then
    Rscript ${script_dir}/figure_setup/sync-salmon-output.R
fi

for project_dir in $salmon_quant_dir/*; do
    project_id=$(basename $project_dir)

    # Step 1: Calculate TPM for each project
    Rscript ${script_dir}/calculate-tpm.R \
        --project_id ${project_id} \
        --output_dir ${tpm_dir}

    # Step 2: Run quanTIseq on each project
    # forthcoming!

done

