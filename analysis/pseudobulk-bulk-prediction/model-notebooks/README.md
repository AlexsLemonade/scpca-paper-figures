This directory contains notebooks which are involved in modeling bulk and pseudobulk.

* `build-assess-models.Rmd` builds and compares fit among sets of candidate linear models for each project
  * It takes a single parameter `expr_threshold` which removes genes which are not expressed in at least that many samples, considering both modalities.
  It exports notebooks named either:
    * `build-assess-models_all-genes.nb.html`: All genes are included (negative threshold specified)
    * `build-assess-models_threshold-{X}.nb.html`: A threshold of `{X}` is specified
