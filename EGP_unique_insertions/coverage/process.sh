#!/bin/bash



#~/workspace/misc_projects/egyptian_refgen/run_egyptref.pl /data/lied_egypt_genome/output_wgs/sample_info.txt /data/lied_egypt_genome/output_wgs2 continue 1 2 1


for f in $(find . -name "*.f13.cov.txt.gz" -exec readlink -f {} \;); do zcat $f | gawk -v v="$f" '{match(v, /\/([^\/]*)\.f13\.cov\.txt\.gz/,a); print v"\t"a[1]"\t"$0}'; done | gzip >>f13.cov.txt.gz


#~/workspace/ngs_analysis/reporting/trimmomatic_report.pl ./ report.trimmomatic.txt
#~/workspace/ngs_analysis/reporting/fastqc_report.pl ./ report.fastqc.txt
#~/workspace/ngs_analysis/reporting/picard_alignment_report.pl ./ report.picard_alignment.txt
