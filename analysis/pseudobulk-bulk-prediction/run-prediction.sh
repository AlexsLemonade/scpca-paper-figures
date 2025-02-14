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
notebook_dir="notebooks"
model_html_dir="${notebook_dir}/model-htmls"
gsea_html_dir="${notebook_dir}/gsea-htmls"

mkdir -p $scpca_dir
mkdir -p $tpm_dir
mkdir -p $pseudobulk_dir
mkdir -p $result_dir
mkdir -p $model_html_dir
mkdir -p $gsea_html_dir

map_file="${data_dir}/bulk-library-sample-ids.tsv"


# Sync data files from S3
Rscript ${script_dir}/sync-data-files.R \
  --output_dir "${scpca_dir}" \
  --map_file "${map_file}"

# Prepare bulk counts data for comparisons
Rscript ${script_dir}/prepare-bulk-counts.R \
  --input_dir "${scpca_dir}" \
  --map_file "${map_file}" \
  --output_counts_file "${data_dir}/normalized-bulk-counts.rds" \
  --output_frac_expressed_file "${data_dir}/fraction-expressed-bulk.tsv"

for project_dir in $scpca_dir/*; do
    project_id=$(basename $project_dir)

    pseudobulk_file="${pseudobulk_dir}/${project_id}_pseudobulk.tsv"
    fraction_expressed_file="${data_dir}/${project_id}_fraction-expressed-single-cell.tsv"

    ###### TPMs are not currently used in the analysis ######
    # Calculate bulk TPM for each project
    #tpm_file="${tpm_dir}/${project_id}_tpm.tsv"
    #Rscript ${script_dir}/calculate-tpm.R \
    #  --input_dir "${scpca_dir}/${project_id}" \
    #  --output_pseudobulk_file "${tpm_file}"

    # Calculate pseudobulk matrices for each project
    Rscript ${script_dir}/calculate-pseudobulk.R \
      --input_dir "${scpca_dir}/${project_id}" \
      --output_pseudobulk_file "${pseudobulk_file}" \
      --output_frac_expressed_file "${fraction_expressed_file}"
done

# Build and export models to results/models across different thresholds for expression
for expr_threshold in -1 0 0.25; do

  if [[ ${expr_threshold} == -1 ]]; then
    threshold_str="all-genes"
  else
    threshold_str="threshold-${expr_threshold}"
  fi

  Rscript -e "rmarkdown::render('${notebook_dir}/build-assess-models.Rmd',
              params = list(expr_threshold = ${expr_threshold}),
              output_file = 'build-assess-models_${threshold_str}.nb.html',
              output_dir = '${model_html_dir}')"
done

# Run the GSEA analysis across gene sets and models
reps=50
for geneset in "H" "C8"; do
  for expr_threshold in -1 0 0.25; do

    if [[ ${expr_threshold} == -1 ]]; then
      threshold_str="all-genes"
    else
      threshold_str="threshold-${expr_threshold}"
    fi

    Rscript -e "rmarkdown::render('${notebook_dir}/perform-gsea.Rmd',
                params = list(msigdbr_category = '$geneset', reps = $reps, model_expr_threshold = ${expr_threshold}),
                output_file = 'perform-gsea_${geneset}_${threshold_str}.nb.html',
                output_dir = '${gsea_html_dir}')"
  done
done
