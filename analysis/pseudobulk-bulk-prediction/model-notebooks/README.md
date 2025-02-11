This directory contains notebooks which are involved in modeling bulk and pseudobulk.

* `build-assess-models.Rmd` builds and compares fit among sets of candidate linear models for each project
  * When run with `../run-prediction.sh`, it produces two notebooks: `build-assess-models_filtered.nb.html` and `build-assess-models_unfiltered.nb.html`.
  In the former, genes not expressed in either modality are removed before modeling each project, and in the latter all genes are retained.
