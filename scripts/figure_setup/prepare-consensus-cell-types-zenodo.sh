#!/bin/bash

# This script creates a tarball `./zenodo-cell-type-consensus-results.tar.gz` with consensus cell type files to upload to Zenodo.
# This script will overwrite existing output, if present.
#
# To run this script, you first need to run `sync-consensus-celltype-results.R`.
# This script is intended for use by the Data Lab only.
#
# Usage: ./prepare-consensus-cell-types-zenodo.sh


# enviroment settings
set -euo pipefail

# Find and go to base directory, which is where this script lives
BASEDIR=$(dirname "${BASH_SOURCE[0]}")
cd "$BASEDIR"

# Define input and output paths
celltype_dir="../../s3_files/cell-type-consensus-results"
output_dir="./zenodo-cell-type-consensus-results"
output_tar="$output_dir.tar.gz"

# Remove existing output directory & tarball, if present
rm -rf $output_dir $output_tar

# Copy the full cell type directory over
cp -r $celltype_dir $output_dir

# Loop over the celltype TSV files to retain only columns needed to reproduce figures and analysis
for celltype_file in ${output_dir}/SCPCP*/SCPCS*/*_processed_consensus-cell-types.tsv.gz; do

  temp_file=$(mktemp) # temp file the tsv.gz gets uncompressed to
  tsv="${celltype_file%.gz}" # file path without .gz for re-compressing at the end

  # Uncompress the tsv file
  gunzip -cf $celltype_file > $temp_file

  # Only retain columns of interest
  cut $temp_file -d$'\t' -f 1,2,3,4,5,12 > $tsv

  # Recompress the tsv file
  gzip -f ${tsv} > ${celltype_file}

  # Remove the temp file
  rm $temp_file
done

# Compress the output directory
tar -czvf $output_tar $output_dir

# Remove the directory
rm -rf $output_dir

