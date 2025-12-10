#!/bin/bash
#
# This bash script runs all scripts to generate manuscript figures and tables.
# This script assumes that you already have locally prepared the necessary input data.
# For full instructions on preparing the input data, please refer to the repository's `README.md` file.
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

# Check if s3 files directory already exists and error if not
if [ ! -d "${BASEDIR}/s3_files" ]; then
  echo "The `s3_files` directory does not exist. To create it, please follow the 'Instructions to prepare data' provided in the README file."
  exit 1
fi

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
Rscript ${script_dir}/Fig2B-F_qc-plots.R

# TODO: this figure will be moving:
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/315
# Figure 4B
Rscript ${script_dir}/Fig4B_singler-cellassign-heatmap.R

# Figure 5A-B
Rscript ${script_dir}/Fig5A-B_nb-umap-panels.R.R

# TODO: this figure will be moving:
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/315
# Figure 5A
# also supplemental figure 6
Rscript ${script_dir}/Fig5A-S6_consensus-cell-type-dotplots.R

# TODO: this figure will be moving:
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/315
# Figures 5B-C
Rscript ${script_dir}/Fig5B-C_brain-barplots.R

# TODO: this figure will be moving:
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/315
# Figure 6
# also supplemental figure 8
Rscript ${script_dir}/Fig6-FigS8_bulk-analysis.R


### Supplementary text figures ###

# Figure S1A
Rscript ${script_dir}/FigS1A_memory-time-comparison.R

# Figure S1B-D
Rscript ${script_dir}/FigS1B-D_method-metrics-comparison.R

# Figure S2B
Rscript ${script_dir}/FigS2B-D_adt-plots.R

# Figure S3D
Rscript ${script_dir}/FigS3D_merged-umaps.R

# TODO: this figure will be moving TBD:
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/315
# Figure S5B
Rscript ${script_dir}/FigS5B_submitter-celltypes-heatmap.R

# TODO: this figure will be moving:
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/315
# Figure S7
Rscript ${script_dir}/FigS7_consensus-bar-plots.R

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
