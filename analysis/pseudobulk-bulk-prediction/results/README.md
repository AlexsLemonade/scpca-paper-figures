This directory contains analysis results.

* `models/`
  * These directory (locally) holds `.rds` files that are ignored from the repository due to large size but can be regenerated with `../run-prediction.sh`.
  * Files are named `<project id>_bulk-pseudobulk-model.rds` and contain a list of two items:
    * `data`: The project-specific data frame, including Ensembl ids, that was modeled
    * `model`: The fitted model itself built as `bulk ~ pseudobulk + (1|sample_id)` with the `lme4` package
