This directory contains analysis results.

* `models/`
  * This directory (locally) holds `.rds` and `.tsv` files that are ignored from the repository due to large size but can be regenerated with `../run-prediction.sh`
  * Files are named according to their corresponding notebooks in `../notebooks/model-htmls`
    * RDS files contain the fitted model itself built as `bulk ~ pseudobulk + (1|sample_id)` with the `lme4` package
    * TSV files contain the underlying data, including Ensembl ids, that was modeled, as well as model residuals
* All TSV files `<project id>_ORA-odds-ratios.tsv` contain results from overrepresentation analysis for the given project
  * Associated HTMLs associated with this analysis are available in `../notebooks/ora-htmls`
