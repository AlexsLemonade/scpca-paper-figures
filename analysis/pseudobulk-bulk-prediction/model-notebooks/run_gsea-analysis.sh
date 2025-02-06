#!/bin/bash

# Run gsea-analysis.Rmd on each project

for project in SCPCP000001 SCPCP000002 SCPCP000006 SCPCP000009 SCPCP000017; do
    output_html="${project}_gsea-analysis.html"
    Rscript -e "rmarkdown::render('gsea-analysis.Rmd', params = list(project_id = '$project'), output_file = '${output_html}')"
done
