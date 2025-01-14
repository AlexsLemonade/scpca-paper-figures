This directory contains results from bulk deconvolution

* `<project id>-quantiseq.tsv` files contain cell type proportions inferred with `quanTIseq` by `scripts/run-quantiseq.R`
* `<project id>-epic.tsv` files contain cell type proportions inferred with `EPIC` by `scripts/run-epic.R`
  * Inferences from both the `Tref` and `Bref` references are included as well as an indicator for whether the model converged for each sample
* `<project id>-epic-full-object.rds` files contain a list of the full `EPIC` objects, named by the reference (`Tref` and `Bref`), inferred by `scripts/run-epic.R`
