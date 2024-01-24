#!/bin/bash
#
# This bash script runs all scripts to generate manuscript figures and tables.
# This script assumes that you already have locally downloaded the necessary data from S3.
# Currently, this requires that you have access to the relevant S3 bucket as a Data Lab member.
#
# To obtain all data, first run the following scripts (paths are relative to the base directory)
#
# Rscript scripts/figure_setup/sync-data-files.R
# Rscript scripts/figure_setup/sync-metadata.R
# Rscript scripts/figure_setup/sync-reference-files.R
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
mkdir -p ${BASEDIR}/tables

#########################################################
#        Generate figures in order of appearance        #
#########################################################


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

# Figure S4A-B
Rscript ${script_dir}/FigS4A-B_celltype-diagnostic-plots.R

##########################################################
#         Generate tables in order of appearance         #
##########################################################


### Supplementary text tables ###

Rscript ${script_dir}/TableS1_modality-summary.R






##########################################################
#                        Clean up                        #
##########################################################

rm -f Rplots.pdf
