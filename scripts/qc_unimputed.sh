#!/bin/bash


files=("cases" "controls")
hwe=( ["cases"]=0.00001 ["controls"]=0.001)
reference="human_g1k_v37.fasta"
ref_af="/matthias/2018-05-28_-_grs_bahnsen/af.vcf.gz"

for f in "${files[@]}"
do

    # Check sex
    plink1.9 --bfile $f --check-sex --out $f


    # Keep autosomes only
    plink --bfile $f --autosome --make-bed --out $f.auto


    # Keep SNPs only and remove AT and GC genotypes
    awk '{if(!(($5 == "A" && $6 == "G") || ($5 == "G" && $6 == "A") || ($5 == "A" && $6 == "C") || ($5 == "C" && $6 == "A") || ($5 == "T" && $6 == "G") || ($5 == "G" && $6 == "T") || ($5 == "T" && $6 == "C") || ($5 == "C" && $6 == "T")))print $0}' $f.auto.bim >$f.snps.remove.txt
    plink --bfile $f.auto --exclude $f.snps.remove.txt --make-bed --out $f.auto.rmsnp


    # General filtering
    plink --bfile $f.auto.rmsnp --geno 0.02 --make-bed --out $f.auto.rmsnp.geno
    plink --bfile $f.auto.rmsnp.geno --hwe ${hwe[$f]}  --make-bed --out $f.auto.rmsnp.geno.hwe
    plink --bfile $f.auto.rmsnp.geno.hwe --mac 10 --make-bed --out $f.auto.rmsnp.geno.hwe.mac


    # Convert to VCF to use bcftools
    plink --bfile $f.auto.rmsnp.geno.hwe.mac --recode vcf id-delim=":" -out $f.auto.rmsnp.geno.hwe.mac


    # Swap and flip alleles according to reference
    bcftools +fixref $f.auto.rmsnp.geno.hwe.mac.vcf -Oz -o $f.auto.rmsnp.geno.hwe.mac.fixref.vcf.gz -- -f $reference -m flip
    bcftools norm --check-ref ws -f $reference $f.auto.rmsnp.geno.hwe.mac.fixref.vcf.gz -Oz -o $f.auto.rmsnp.geno.hwe.mac.fixref.norm.vcf.gz


    # Sort, recalculate stats, create index
    bcftools sort -Oz $f.auto.rmsnp.geno.hwe.mac.fixref.norm.vcf.gz $f.auto.rmsnp.geno.hwe.mac.fixref.norm.vcf.gz | bcftools +fill-tags -Oz -o $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.vcf.gz 
    bcftools index $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.vcf.gz


    # Check genotype frequency and check high alt allele frequencies
    bcftools annotate -c INFO/AF -a $ref_af $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.vcf.gz | bcftools +af-dist | grep ^PROB >$f.af-dist.txt
	~/workspace/gwas/pipeline/functions/vis_afdist.pl $f.af-dist.txt $f.af-dist.png
    bcftools view -i'INFO/AF>0.95' $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.vcf.gz -Ov -o $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.hialtfreq.vcf.gz


    # Convert back to plink
    plink --vcf $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.vcf.gz --make-bed --out $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann
    cut -f1,2,5 $f.auto.rmsnp.fam >$f.sex.txt
    cut -f1,2,6 $f.auto.rmsnp.fam >$f.pheno.txt
    plink1.9 --bfile $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann --update-sex $f.sex.txt --make-pheno $f.pheno.txt '2' --make-bed --out $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd --allow-no-sex --keep-allele-order


    # Normalize SNP identifier
    ~/workspace/plink/normalize_snp_ids_in_plink.pl $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.bim $f.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.norm_ids.bim $f.id_mappings.txt 100 0
done


# Merge cases and controls
plink1.9 --fam cases.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.fam \
        --bed cases.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.bed \
        --bim cases.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.norm_ids.bim \
        --bmerge controls.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.bed \
               controls.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.norm_ids.bim \
                 controls.auto.rmsnp.geno.hwe.mac.fixref.norm.ann.upd.fam --make-bed --out all --keep-allele-order --allow-no-sex


# Keep SNPs that are available in all datasets
plink --bfile all --geno 0.02 --make-bed --out all.geno


# Remove individuals with low genotype call rate
plink --bfile all.geno --mind 0.95 --make-bed --out all.geno.mind


# Remove by het rate
plink1.9 --bfile all.geno.mind --het --out all.geno.mind --keep-allele-order --allow-no-sex
~/workspace/gwas/pipeline/functions/find_het_outliers.pl all.geno.mind.het hetfail.txt 3
awk '{if($6>0.2 || $6<-0.2) print $1 " " $2}' all.geno.mind.het >>hetfail.txt
plink --bfile all.geno.mind --remove hetfail.txt --make-bed --out all.geno.mind.hetfail


# Find cryptic related individuals (LD pruning is not recommended in KING)
# >0.0884 => remove 2nd degree and closer
plink --bfile all.geno.mind.hetfail --make-king-table --king-cutoff 0.0884 --make-bed --out all.geno.mind.hetfail.king


# Find outliers by PCA
/usr/local/bin/python3.7 ~/workspace/gwas/pipeline/functions/pc_outlier_rm.py --npc-rm 5 --nsd-rm 5 --bfile all.geno.mind.hetfail.king --out all.geno.mind.hetfail.king.pca
cat ~/workspace/plink/merge/plot_pc.R | R --slave --args all.geno.mind.hetfail.king.fam all.geno.mind.hetfail.king.eigenvec pca.scatterplot.before.png
cat ~/workspace/plink/merge/plot_pc.R | R --slave --args all.geno.mind.hetfail.king.pca.fam all.geno.mind.hetfail.king.pca.eigenvec pca.scatterplot.after.png


# Filter by MAF
plink --bfile all.geno.mind.hetfail.king.pca --maf 0.05 --make-bed --out all.geno.mind.hetfail.king.pca.maf


# Perform association test
plink --bfile all.geno.mind.hetfail.king.pca --extract all.geno.mind.hetfail.king.prune.in --pca 5 --out all.geno.mind.hetfail.king.pca
plink --bfile all.geno.mind.hetfail.king.pca.maf --glm sex genotypic --out all.geno.mind.hetfail.king.pca.maf --ci .95 --covar-name PC1-PC5 --covar all.geno.mind.hetfail.king.pca.eigenvec


Manhattan plot
Rscript -e 'library(data.table); dat = fread("all.geno.mind.hetfail.king.pca.maf.PHENO1.glm.logistic"); chisq=qchisq(1-dat$P, df=1); lambda=median(chisq)/qchisq(0.5,1); print(paste0("Lambda: ", lambda));'
Rscript -e 'library(qqman); library(data.table); dat = fread("all.geno.mind.hetfail.king.pca.maf.PHENO1.glm.logistic"); names(dat)[1:2]=c("CHR", "BP"); png(filename="glm.manhattan2.png", units="mm", res=1200, width=297, height=210, bg="transparent"); manhattan(dat[!is.na(P) & TEST == "ADD"]); dev.off()'
Rscript -e 'library(qqman); library(data.table); dat = fread("all.geno.mind.hetfail.king.pca.maf.PHENO1.glm.logistic"); names(dat)[1:2]=c("CHR", "BP"); png(filename="glm.qq.png", units="mm", res=1200, width=297, height=210, bg="transparent"); qq(dat[!is.na(P) & TEST == "ADD"]$P); dev.off()'


plink --bfile all.geno.mind.hetfail.king.pca --chr 6 -make-bed --out all.geno.mind.hetfail.king.pca.chr6
awk 'BEGIN{print "chr pos snp ea or l95 u95 p"} NR > 1 {if($7 == "ADD" && $14 != "NA") print $1 " " $2 " " $3 " " $6 " " $9 " " $11 " " $12 " " $14}' all.geno.mind.hetfail.king.pca.maf.PHENO1.glm.logistic >assoc_results.txt
~/workspace/snp_list_tools/cluster.pl assoc_results.txt 200000 0.0001 1 chr pos p assoc_results.cluster.p.txt assoc_results.cluster.pos.txt assoc_results.cluster.txt


# Imputation ready
plink1.9 --make-bed --bfile all.geno.mind.hetfail.king.pca --filter-cases -out all.geno.mind.hetfail.king.pca.cases --keep-allele-order --allow-no-sex
plink1.9 --make-bed --bfile all.geno.mind.hetfail.king.pca --filter-controls -out all.geno.mind.hetfail.king.pca.controls --keep-allele-order --allow-no-sex
plink --bfile all.geno.mind.hetfail.king.pca.cases --recode vcf id-delim=":" --out all.geno.mind.hetfail.king.pca.cases
plink --bfile all.geno.mind.hetfail.king.pca.controls --recode vcf id-delim=":" --out all.geno.mind.hetfail.king.pca.controls
bcftools +fixref all.geno.mind.hetfail.king.pca.cases.vcf -Oz -o all.geno.mind.hetfail.king.pca.cases.vcf.gz -- -f $reference -m flip
bcftools +fixref all.geno.mind.hetfail.king.pca.controls.vcf -Oz -o all.geno.mind.hetfail.king.pca.controls.vcf.gz -- -f $reference -m flip


# Check rare variants for MAF errors
plink1.9 --bfile all.geno.mind.hetfail.king.pca  --assoc --out all.geno.mind.hetfail.king.pca --allow-no-sex
awk '{if($9 == "NA" || $9 < 1E-10) print $0}' all.geno.mind.hetfail.king.pca.assoc >snps.rm_assoc.txt


# SNPtest
awk '{print $2 " " $2 " " $6-1}' all.geno.mind.hetfail.king.pca.fam >pheno.txt
awk '{$1=$2; print$0}' all.geno.mind.hetfail.king.pca.eigenvec | cut -f1-7 -d " " >cov.eigenvec.txt
awk '{print $2 " " $2 " " $5}' all.geno.mind.hetfail.king.pca.fam >sex.txt
