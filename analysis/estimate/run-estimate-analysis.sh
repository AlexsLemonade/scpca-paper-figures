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
tpm_dir="${data_dir}/tpm"
pseudobulk_dir="${data_dir}/pseudobulk"

mkdir -p $scpca_dir
mkdir -p $tpm_dir
mkdir -p $pseudobulk_dir

library_sample_map_file="${data_dir}/bulk-library-sample-ids.tsv"
ensembl_symbol_map_file="${data_dir}/ensembl-symbol-map.tsv"

# Sync data files from S3
Rscript ${script_dir}/sync-data-files.R \
  --output_dir "${scpca_dir}" \
  --library_sample_map_file "${library_sample_map_file}" \
  --ensembl_symbol_map_file "${ensembl_symbol_map_file}"
