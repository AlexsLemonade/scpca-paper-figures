#!/bin/bash
#
# This bash script runs all scripts to generate manuscript figures and tables.
# This script assumes that you already have locally downloaded the necessary data from S3.
# Currently, this requires that you have access to the relevant S3 bucket as a Data Lab member.
#
# This script can be run as:
# bash generate-figures-tables.sh

# enviroment settings
set -euo pipefail

# Find and go to base directory, which is where this script lives
BASEDIR=$(dirname "${BASH_SOURCE[0]}")
cd "$BASEDIR"

# Define path to figure/table generation scripts
script_dir=${BASEDIR}/scripts

# Ensure output directories exist
mkdir -p ${BASEDIR}/figures/pngs
mkdir -p ${BASEDIR}/figures/pdfs
mkdir -p ${BASEDIR}/tables

# check if s3 files directory already exists
# if it does, delete it
if [ -d "${BASEDIR}/s3_files" ]; then
  rm -r ${BASEDIR}/s3_files
fi
mkdir -p ${BASEDIR}/s3_files

# copy over files from S3
echo "Downloading data files..."
Rscript ${script_dir}/figure_setup/sync-data-files.R
Rscript ${script_dir}/figure_setup/sync-metadata.R
Rscript ${script_dir}/figure_setup/sync-reference-files.R

#########################################################
#        Generate figures in order of appearance        #
#########################################################

echo "Generating figures and tables..."


### Main text figures ###

# Figure 1A
Rscript ${script_dir}/Fig1A_sample-disease-barchart.R

# Figure 1B
Rscript ${script_dir}/Fig1B_modality-barchart.R

# Figure 2B
Rscript ${script_dir}/Fig2B_qc-plots.R

# Figure 4B
Rscript ${script_dir}/Fig4B_singler-cellassign-heatmap.R

### Supplementary text figures ###

# Figure S1A
Rscript ${script_dir}/FigS1A_memory-time-comparison.R

# Figure S1B-D
Rscript ${script_dir}/FigS1B-D_method-metrics-comparison.R

# Figure S2B
Rscript ${script_dir}/FigS2B_adt-plots.R

# Figure 3D
Rscript ${script_dir}/Fig3D_merged-umaps.R

# Figure S4A-B
Rscript ${script_dir}/FigS4A-B_celltype-diagnostic-plots.R

# Figure S5
Rscript ${script_dir}/FigS5_submitter-celltypes-heatmap.R

##########################################################
#         Generate tables in order of appearance         #
##########################################################


### Supplementary text tables ###

Rscript ${script_dir}/TableS1_modality-summary.R

Rscript ${script_dir}/TableS2_cellassign-ref-summary.R




##########################################################
#                        Clean up                        #
##########################################################

rm -f Rplots.pdf
