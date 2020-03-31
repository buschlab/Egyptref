# GWAS
This repository provides a description of workflow for analysing GWAS data. All the required scripts are also included. A GWAS analysis can be subdivided into 5 major steps:

1. QC of unimputed data
2. Imputation
3. QC of imputed data
4. Association testing
5. Regional clustering and plotting

Note: In some cases, scripts, parameters and thresholds need to be adapted to consider certain data proporties.

## 1. QC of unimputed data
**Script:** qc_unimputed.sh \
**Input data:** Seperate Plink files for cases and controls with sex info and both family id and individual id set \
**Required programs:** plink2.0, plink1.9, bcftools, python3.7, Rscript with libraries qqman and data.table \
**Annotation files:** reference genome (fasta)


## 2. Imputation
For imputation I recommend the Sanger imputation server (free of charge) at https://imputation.sanger.ac.uk/. Either select "pre-phase with EAGLE2 and impute" or "pre-phase with SHAPEIT2 and impute" and upload the final VCF files from step 1.


## 3. QC of imputed data
