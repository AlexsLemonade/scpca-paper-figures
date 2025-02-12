This directory contains exploratory notebook that are not directly used in the analysis pipeline.

* `compare-pseudobulk-bulk.Rmd`: Compares different approaches of calculating pseudobulk and bulk expression
* `residuals-annotated-genes.Rmd`: Explores distribution of model residuals between genes with and without official gene symbols
* `permute-gsea.Rmd`: Performs GSEA on model residuals, using shuffled version of a given `MSigDB` signature set over a given number of replicates
  * This serves as a null test of our use of this approach to identify differences between modalities
  * When run with `../run_permute-gsea.sh`, it outputs these notebooks to `permute-gsea-html`:
    * `permute-gsea_{H, C8}_{all-genes, threshold-0, threshold-0.25}.nb.html`
    * The first `{field}` refers to the `MSigDB` signature set, and the second `{field}`  refers to which model's residuals the notebook used; see `../model-notebooks/README.md` for more context

