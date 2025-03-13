This directory contains scripts and results for predicting bulk expression from matched pseudobulk expression.

To run the analysis, run:

```sh
bash run-prediction.sh
```

If you are a Data Lab member using 1Password to manage your credentials, you may need to run:
```sh
op run -- bash run-prediction.sh
```

The analysis takes the following steps:

* Transform bulk count matrices with `DESeq2`
* Prepare pseudobulk matrices and transform with `DESeq2`
* Prepare consensus cell types gene sets for use in overrepresentation analysis
* Construct linear models predicting bulk from pseudobulk expression
* Perform overrepresentation analysis on residual outliers
