#!/bin/bash


path='/data/lied_egypt_genome/output_wgs/sv'



cat <(echo $'chr\tpos\tfilter\tend\ttsv_id\tsv_type\tsample\tgt\tstatus') \
	<(bcftools query -i 'FILTER="PASS"' -f '[%CHROM\t%POS\t%FILTER\t%INFO/END\t%ID\t%INFO/SVTYPE\t%SAMPLE\t%GT\t%FT\n]' egyptians_deletions.vcf \
	| grep -v "0/0" | grep PASS) >del.txt


cat <(echo $'chr\tpos\tfilter\tend\ttsv_id\tsv_type\tsample\tgt\tstatus') \
        <(bcftools query -i 'FILTER="PASS"' -f '[%CHROM\t%POS\t%FILTER\t%INFO/END\t%ID\t%INFO/SVTYPE\t%SAMPLE\t%GT\t%FT\n]' egyptians_inversions.vcf \
        | grep -v "0/0" | grep PASS) >inv.txt


cat <(echo $'chr\tpos\tfilter\tend\ttsv_id\tsv_type\tsample\tgt\tstatus') \
        <(bcftools query -i 'FILTER="PASS"' -f '[%CHROM\t%POS\t%FILTER\t%INFO/END\t%ID\t%INFO/SVTYPE\t%SAMPLE\t%GT\t%FT\n]' egyptians_duplications.vcf \
        | grep -v "0/0" | grep PASS) >dupl.txt


cat <(echo $'chr\tpos\tfilter\tend\ttsv_id\tsv_type\tsample\tgt\tstatus\tinsertion_length') \
        <(bcftools query -i 'FILTER="PASS"' -f '[%CHROM\t%POS\t%FILTER\t%INFO/END\t%ID\t%INFO/SVTYPE\t%SAMPLE\t%GT\t%FT\t%INFO/INSLEN\n]' egyptians_insertions.vcf \
        | grep -v "0/0" | grep PASS) >ins.txt


cat <(echo $'chr\tpos\tfilter\tchr2\tend\ttsv_id\tsv_type\tsample\tgt\tstatus') \
	<(bcftools query -i 'FILTER="PASS"' -f '[%CHROM\t%POS\t%FILTER\t%INFO/CHR2\t%INFO/END\t%ID\t%INFO/SVTYPE\t%SAMPLE\t%GT\t%FT\n]' egyptians_translocations.vcf \
        | grep -v "0/0" | grep PASS) >transloc.txt


#https://stackoverflow.com/questions/16957293/collapse-intersecting-regions
#https://groups.google.com/forum/#!msg/delly-users/wD1xZXk8Jt8/eH1SZiyj1MoJF

