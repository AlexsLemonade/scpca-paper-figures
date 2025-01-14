This directory contains results from bulk deconvolution

* `<project id>-quantiseq.tsv` files contain cell type proportions inferred with `quanTIseq` by `scripts/run-quantiseq.R`
* `<project id>-epic.tsv` files contain cell type proportions inferred with `EPIC` by `scripts/run-epic.R`
  * Includes inferences from both the `TRef` and `BRef` references
  * Includes a column `mRNAProportions` which contains cell type fractions _without_ correcting for mRNA content among cell types
  * Includes an indicator column for whether the model converged for each sample
