# GWAS
This repository provides a description of workflow for analysing GWAS data. All the required scripts are also included. A GWAS analysis can be subdivided into 5 major steps:

1. QC of unimputed data
2. Imputation
3. QC of imputed data
4. Association testing
5. Regional clustering and plotting


## 1. QC of unimputed data
**Script:** qc_unimputed.sh \
**Input data:** Seperate Plink files for cases and controls with sex info and both family id and individual id set \
**Required programs:** plink2.0, plink1.9, bcftools, python3.7, Rscript with libraries qqman and data.table \
**Annotation files:** reference genome (fasta)


## 2. Imputation
