#!/bin/bash

path="..."

#~/workspace/gwas/snptest/excl_all_position_duplicates.pl dups.txt ../cases_1000GP3/gen/*.gen.gz ../controls_1000GP3/gen/*.gen.gz

#cat <(head -n2 cov.eigenvec.txt) <(sed 1,2d cov.eigenvec.txt | awk '{id=$1 ":" $2; $1=id; $2=id; print $0}') >cov.txt
#cat <(head -n2 all.geno.maf.fam | cut -f1,2,6) <(sed 1,2d all.geno.maf.fam | awk '{id=$1 ":" $2; $1=id; $2=id; print $1 "\t" $2 "\t" $6}') >pheno.txt
#~/workspace/gen/add_covariates.pl cov.eigenvec.txt 1.samples cases.samples

~/workspace/gwas/snptest/run_snptest.pl prefixes.txt dups.txt pheno "sex PC1 PC2 PC3 PC4 PC5" 22

mkdir results
~/workspace/toolbox/fork.pl 11 \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr1.out $path/results/snptest.chr1.txt add 1 0. 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr2.out $path/results/snptest.chr2.txt add 2 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr3.out $path/results/snptest.chr3.txt add 3 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr4.out $path/results/snptest.chr4.txt add 4 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr5.out $path/results/snptest.chr5.txt add 5 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr6.out $path/results/snptest.chr6.txt add 6 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr7.out $path/results/snptest.chr7.txt add 7 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr8.out $path/results/snptest.chr8.txt add 8 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr9.out $path/results/snptest.chr9.txt add 9 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr10.out $path/results/snptest.chr10.txt add 10 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr11.out $path/results/snptest.chr11.txt add 11 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr12.out $path/results/snptest.chr12.txt add 12 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr13.out $path/results/snptest.chr13.txt add 13 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr14.out $path/results/snptest.chr14.txt add 14 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr15.out $path/results/snptest.chr15.txt add 15 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr16.out $path/results/snptest.chr16.txt add 16 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr17.out $path/results/snptest.chr17.txt add 17 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr18.out $path/results/snptest.chr18.txt add 18 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr19.out $path/results/snptest.chr19.txt add 19 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr20.out $path/results/snptest.chr20.txt add 20 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr21.out $path/results/snptest.chr21.txt add 21 0 0.8" \
"~/workspace/gwas/snptest/create_result_file.pl $path/chr22.out $path/results/snptest.chr22.txt add 22 0 0.8"


~/workspace/toolbox/file_tools/concat.pl $path/results/snptest.all.txt 1 $path/results/snptest.chr*.txt
