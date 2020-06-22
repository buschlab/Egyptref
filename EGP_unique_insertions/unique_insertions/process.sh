###!/bin/bash


path='/data/lied_egypt_genome/output_wgs2/novel_sequences'

#module load Perl-Packages/5
#module load URI/1.76
#module load LWP-UserAgent/6.39
#module load Statistics-Basic/1.6611

#for f in $(find ../ -name "*.f13.cov.txt.gz" -exec readlink -f {} \;); do zcat $f | gawk -v v="$f" '{match(v, /\/([^\/]*)\.f13\.cov\.txt\.gz/,a); print v"\t"a[1]"\t"$0}'; done | gzip >>f13.cov.txt.gz
#~/workspace/misc_projects/egyptian_refgen/transform_cov.pl f13.cov.txt f13.cov.transformation.txt

#~/workspace/misc_projects/egyptian_refgen/extract_high_cov_seqs.pl f13.cov.transformation.txt novel_seq.gr5_50 \
#/data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa n_gr_5 50 500

#~/workspace/misc_projects/egyptian_refgen/extract_high_cov_seqs.pl f13.cov.transformation.txt novel_seq.gr20_10 \
#/data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa n_gr_20 10 500

#~/workspace/misc_projects/egyptian_refgen/extract_high_cov_seqs.pl f13.cov.transformation.txt novel_seq.gr5_10 \
#/data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa n_gr_5 10 500

#~/workspace/misc_projects/egyptian_refgen/extract_high_cov_seqs.pl f13.cov.transformation.txt novel_seq.gr1_1 \
#/data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa n_gr_zero 1 500


## Variant calling
#for f in $(find ../ -name "*.f13.bam" -exec readlink -f {} \;); do
#	echo $f >>bam_files.txt
#done

#grep -h ">" seqs/novel_seq.gr5_50.*.fasta | perl -n -e '/>(.*):([0-9]+)-([0-9]+)$/ && print $1."\t".$2."\t".$3."\n"' >novel_seq.gr5_50.regions.txt
#grep -h ">" seqs/novel_seq.gr20_10.*.fasta | perl -n -e '/>(.*):([0-9]+)-([0-9]+)$/ && print $1."\t".$2."\t".$3."\n"' >novel_seq.gr20_10.regions.txt
#grep -h ">" seqs/novel_seq.gr5_10.*.fasta | perl -n -e '/>(.*):([0-9]+)-([0-9]+)$/ && print $1."\t".$2."\t".$3."\n"' >novel_seq.gr5_10.regions.txt
#grep -h ">" seqs/novel_seq.gr1_1.*.fasta | perl -n -e '/>(.*):([0-9]+)-([0-9]+)$/ && print $1."\t".$2."\t".$3."\n"' >novel_seq.gr1_1.regions.txt


#samtools mpileup -u -g -l novel_seq.gr5_50.regions.txt -b bam_files.txt -f /data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa \
#| bcftools call -v -m -O z -o variants/novel_seq.gr5_50.vcf.gz

#samtools mpileup -u -g -l novel_seq.gr20_10.regions.txt -b bam_files.txt -f /data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa \
#| bcftools call -v -m -O z -o variants/novel_seq.gr20_10.vcf.gz

#samtools mpileup -u -g -l novel_seq.gr5_10.regions.txt -b bam_files.txt -f /data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa \
#| bcftools call -v -m -O z -o variants/novel_seq.gr5_10.vcf.gz

#samtools mpileup -u -g -l novel_seq.gr1_1.regions.txt -b bam_files.txt -f /data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa \
#| bcftools call -v -m -O z -o variants/novel_seq.gr1_1.vcf.gz


for f in $(ls variants/*vcf.gz); do
	bcftools stats -F /data/lied_egypt_genome/reference/EGYPTREFMETAV2ADDED/Homo_sapiens.EGYPTREFMETAV2ADDED.dna.primary_assembly.fa $f >$f.stats
	mkdir ${f}_plots
	plot-vcfstats -p ${f}_plots/ $f.stats
done
