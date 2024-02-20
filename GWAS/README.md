**README for data_wrangle.R & GWAS.R**

**data_wrangle.R**
-Load neuronal and nonneuronal datasets
-Normalized the raw counts table to be 10,000 before calculate quantile groups
-Calculate the quantile groups for each celltype


**GWAS.R**
-Load quantile groups for celltypes saved from data_wrangle.R
-Check whether gwas is properly formatted
-Map SNPs to Genes
-Test GWAS enrichments that remain in a GWAS after controling for a second GWAS


