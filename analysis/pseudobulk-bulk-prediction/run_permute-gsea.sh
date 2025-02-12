#!/bin/bash

# This script runs the GSEA shuffling notebook under several conditions

set -euo pipefail

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd "$basedir"

# Define directories
notebook_dir="exploratory-notebooks"

reps=50 # 100 would be nice but we have places to be

# Run across gene sets and filtered versions
for geneset in "H" "C8"; do
  for expr_threshold in -1 0 0.25; do

    if [[ ${expr_threshold} == -1 ]]; then
      threshold_str="all-genes"
    else
      threshold_str="threshold-${expr_threshold}"
    fi

    Rscript -e "rmarkdown::render('${notebook_dir}/permute-gsea.Rmd',
                params = list(msigdbr_category = '$geneset', reps = $reps, model_expr_threshold = ${expr_threshold}),
                output_file = 'permute-gsea_${geneset}_${threshold_str}.nb.html',
                output_dir = '${notebook_dir}/permute-gsea-html')"
  done
done

