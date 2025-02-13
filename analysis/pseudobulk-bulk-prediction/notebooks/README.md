This directory contains notebooks which are used for analysis

* `build-assess-models.Rmd` builds and compares fit among sets of candidate linear models for each project
  * It takes a single parameter `expr_threshold` which removes genes which are not expressed in at least that many samples, considering both modalities.
  It exports notebooks named either:
    * `build-assess-models_all-genes.nb.html`: All genes are included (negative threshold specified)
    * `build-assess-models_threshold-{X}.nb.html`: A threshold of `{X}` is specified
* `perform-gsea.Rmd` runs GSEA using an `MSigDB` gene signature set on model residuals,
  * It additionally performs GSEA using a shuffled version of the same gene set to serve as a form of null test of our use of this approach
  * It exports notebooks named as `permute-gsea_{H, C8}_{all-genes, threshold-0, threshold-0.25}.nb.html` to the directory `gsea-notebooks/`
    * The first `{field}` refers to the `MSigDB` signature set, and the second `{field}`  refers to which model's residuals the notebook used
