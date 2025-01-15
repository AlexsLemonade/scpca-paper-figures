This directory contains notebooks used for exploratory analyses.

* `quantiseq-tumor-genes.Rmd` explores expression distributions of genes in `quanTIseq`'s signature gene set among ScPCA samples to determine whether the `tumor=TRUE` setting should be used on ScPCA data
  * This notebook also identifies discrepancies between our gene symbols and those used in `quanTIseq`'s gene signature set and determines how to recode them for analysis
* `explore-quantiseq-results.Rmd` explores, primarily using visualization, cell type distributions across samples and projects inferred by `quanTIseq`
* `epic-signature-genes.Rmd` identifies discrepancies between our gene symbols and those used in `EPIC`'s gene signature sets and determines how to recode them for analysis
* `explore-epic-results.Rmd` explores, primarily using visualization, cell type distributions across samples and projects inferred by `EPIC`'s two references, and identifies which reference should be used moving forward
