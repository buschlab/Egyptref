#!/bin/bash


~/workspace/gwas/pipeline/after_imputation/01_afdist.pl raw/ afdist 20
~/workspace/gwas/pipeline/after_imputation/01_ic.pl raw/ ic
~/workspace/gwas/pipeline/after_imputation/02_filter_info.pl raw/ filter 22 0.8 0.01
~/workspace/gwas/pipeline/after_imputation/03_vcf2gen.pl filter gen 22
