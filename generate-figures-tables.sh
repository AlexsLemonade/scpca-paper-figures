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

# Figure 3C
Rscript ${script_dir}/Fig3C_celltype-methods-heatmap.R

# Figure 4A and S5
Rscript ${script_dir}/Fig4A-S5_consensus-cell-type-dotplots.R

# Figures 4B-C
# TODO: this script needs to be re-run as described in this issue
# https://github.com/AlexsLemonade/scpca-paper-figures/issues/297
Rscript ${script_dir}/Fig4B-C_brain-barplots.R

# Figure 4D 
Rscript ${script_dir}/Fig4D_immune-only-dotplot.R

# Figure 5A-B
Rscript ${script_dir}/Fig5A-B_nb-umap-panels.R

# Figure 5C
Rscript ${script_dir}/Fig5C_nb-cnv-heatmaps.R 

# Figure 5D
Rscript ${script_dir}/Fig5D_nb-ridgeplot.R 

# Figure 6 and S7
Rscript ${script_dir}/Fig6-FigS7_bulk-analysis.R


### Supplementary text figures ###

# Figure S1A
Rscript ${script_dir}/FigS1A_memory-time-comparison.R

# Figure S1B-D
Rscript ${script_dir}/FigS1B-D_method-metrics-comparison.R

# Figure S2B
Rscript ${script_dir}/FigS2B-D_adt-plots.R

# Figure S3D
Rscript ${script_dir}/FigS3D_merged-umaps.R

# Figure S4
Rscript ${script_dir}/FigS4_celltype-method-umaps.R

# Figure S4
# TODO: this script represents an older version of the heatmap, and it may be updated to use in S4
# Or, it may not be used since we may be changing this figure to UMAPs as described:
#https://github.com/AlexsLemonade/scpca-paper-figures/issues/299
Rscript ${script_dir}/FigS4_submitter-celltypes-heatmap.R

# Figure S6
Rscript ${script_dir}/FigS6_consensus-bar-plots.R

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
