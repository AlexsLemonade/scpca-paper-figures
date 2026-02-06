#!/bin/bash

# This script runs the analysis predicting bulk from pseudobulk followed by overrepresentation analysis on model outliers
#
# Default usage:
# ./run-prediction.sh
#
# To also run the GSEA analysis (not presented in the manuscript), use:
# run_gsea=1 ./run-prediction.sh
# 
# To write model HTMLs considering a different threshold for expressed genes across samples to consider, use:
# expr_threshold=0.25 ./run-prediction.sh # for a 0.25 threshold
# expr_threshold=-1 ./run-prediction.sh   # for no threshold


set -euo pipefail

# controls whether to run the GSEA analysis, which is very time-consuming and not ultimately part of the paper
run_gsea=${run_gsea:-0}

# controls the expression threshold to use when modeling; default (used in paper) is 0
expr_threshold=${expr_threshold:-0}

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd "$basedir"

#### Define directories

# Data directories
data_dir="data"
scpca_dir="${data_dir}/scpca_data"
ref_dir="${data_dir}/references"
pseudobulk_dir="${data_dir}/pseudobulk"
consensus_celltype_dir="../../s3_files/cell-type-consensus-results"

# Code directories
script_dir="scripts"
notebook_dir="notebooks"

# Output directories
result_dir="results"
model_html_dir="${notebook_dir}/model-htmls"
gsea_html_dir="${notebook_dir}/gsea-htmls"
ora_html_dir="${notebook_dir}/ora-htmls"


# Check if input directories already exists and error if not
if [[ ! -d $scpca_dir ]] || [[ ! -d $ref_dir ]] || [[ ! -d $consensus_celltype_dir ]] ; then
  echo "At least one input data directory does not exist. To prepare input data, please follow the 'Instructions to prepare data' provided in the reproduce-figures/README file."
  exit 1
fi

# set up usage of expr_threshold
if [[ ${expr_threshold} == -1 ]]; then
  threshold_str="all-genes"
else
  threshold_str="threshold-${expr_threshold}"
fi

# Create additional directories
mkdir -p $pseudobulk_dir
mkdir -p $result_dir
mkdir -p $model_html_dir
mkdir -p $ora_html_dir

# This convenience file keeps track of the bulk library & sample ids used
# This is also used to assist filtering to single-cell processed RDS files of interest based on sample
map_file="${data_dir}/bulk-library-sample-ids.tsv"

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

    # Calculate pseudobulk matrices for each project
    Rscript ${script_dir}/calculate-pseudobulk.R \
      --input_dir "${scpca_dir}/${project_id}" \
      --map_file "${map_file}" \
      --output_pseudobulk_file "${pseudobulk_file}" \
      --output_frac_expressed_file "${fraction_expressed_file}"

done

# Build and export models to results/models at a 0 threshold
Rscript -e "rmarkdown::render(
  '${notebook_dir}/build-assess-models.Rmd',
  params = list(expr_threshold = ${expr_threshold}),
  output_file = 'build-assess-models_${threshold_str}.nb.html',
  output_dir = '${model_html_dir}'
)"



# If specified, run the GSEA analysis across gene sets and models
if [[ ${run_gsea} -eq 1 ]]; then
  mkdir -p $gsea_html_dir

  gsea_reps=50
  for geneset in "H" "C8"; do
    Rscript -e "rmarkdown::render('${notebook_dir}/perform-gsea.Rmd',
                params = list(msigdbr_category = '$geneset', reps = ${gsea_reps}, model_expr_threshold = ${expr_threshold}),
                output_file = 'perform-gsea_${geneset}_${threshold_str}.nb.html',
                output_dir = '${gsea_html_dir}')"
  done
fi

# Run the overrepresentation analysis across gene sets using the model with genes present in at least one modality per sample
ora_reps=10000 # replicates for permutation p-value calculation
summary_function="median" # use median of residuals when summarizing project
sd_threshold=2.5 # outliers are >2.5 sd
for project_dir in $scpca_dir/*; do
  project_id=$(basename $project_dir)
  
    # Determine Panglao file with gene sets use for ORA
    case ${project_id} in
         "SCPCP000001" | "SCPCP000002" | "SCPCP000009")
          panglao_file="brain-compartment_PanglaoDB_2020-03-27.tsv"
          ;;
        "SCPCP000006")
          panglao_file="kidney-compartment_PanglaoDB_2020-03-27.tsv"
          ;;
        "SCPCP000017")
          panglao_file="bone-and-soft-tissue_PanglaoDB_2020-03-27.tsv"
          ;;
      esac

  Rscript -e "rmarkdown::render('${notebook_dir}/perform-ora.Rmd',
                   params = list(project_id = '$project_id',
                                 reps = ${ora_reps},
                                 summary_function = '${summary_function}',
                                 sd_threshold = ${sd_threshold},
                                 panglao_file = '${panglao_file}'),
                   output_file = 'perform-ora_${project_id}.nb.html',
                   output_dir = '${ora_html_dir}')"
done


# render the quick visualization notebook
Rscript -e "rmarkdown::render('exploratory-notebooks/view-ora-results.Rmd')"
