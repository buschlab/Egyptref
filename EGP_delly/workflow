IND=bamfile_directory/
OUT=genotypes_delly/
REF=Homo_sapiens_assembly38.fasta
BLP=_sbatch_dummy_sv
EXC=human.hg38.excl.tsv

DELLY=delly_v0.8.1

LIST=`ls $IND | grep EGAN`
LIST="$LIST LU18 LU19 LU2 LU22 LU23 LU9 PD114 PD115 PD82 EGYPTREF"

mkdir $OUT

# call SNVs
for i in $LIST; do
	echo RUNNING on $i
	CALL=$(echo $DELLY call --genome $REF $IND$i/$i.merged.mark_dups.base_recal.bam --exclude $EXC --vcffile sites.bcf --outfile $OUT/genotype.$i.bcf )
	echo -e "$(cat $BLP)\n\n$CALL\n" > sbatch_script
	sbatch sbatch_script
	rm sbatch_script
done

# Postprocessing
cd genotypes_delly/
bcftools merge -m id -O b -o ../merged.bcf genotype.EGAN00001101667.bcf genotype.EGAN00001101668.bcf genotype.EGAN00001101669.bcf genotype.EGAN00001101670.bcf genotype.EGAN00001101671.bcf genotype.EGAN00001101672.bcf genotype.EGAN00001101676.bcf genotype.EGAN00001101677.bcf genotype.EGAN00001101678.bcf genotype.EGAN00001101679.bcf genotype.EGAN00001101680.bcf genotype.EGAN00001101681.bcf genotype.EGAN00001101682.bcf genotype.EGAN00001101687.bcf genotype.EGAN00001101688.bcf genotype.EGAN00001101689.bcf genotype.EGAN00001101690.bcf genotype.EGAN00001101692.bcf genotype.EGAN00001101694.bcf genotype.EGAN00001101699.bcf genotype.EGAN00001101700.bcf genotype.EGAN00001101702.bcf genotype.EGAN00001101705.bcf genotype.EGAN00001101706.bcf genotype.EGAN00001101711.bcf genotype.EGAN00001101712.bcf genotype.EGAN00001101713.bcf genotype.EGAN00001101716.bcf genotype.EGAN00001101717.bcf genotype.EGAN00001101718.bcf genotype.EGAN00001101719.bcf genotype.EGAN00001101723.bcf genotype.EGAN00001101724.bcf genotype.EGAN00001101725.bcf genotype.EGAN00001101732.bcf genotype.EGAN00001101734.bcf genotype.EGAN00001101735.bcf genotype.EGAN00001101736.bcf genotype.EGAN00001101737.bcf genotype.EGAN00001101739.bcf genotype.EGAN00001101742.bcf genotype.EGAN00001101744.bcf genotype.EGAN00001101748.bcf genotype.EGAN00001101749.bcf genotype.EGAN00001101750.bcf genotype.EGAN00001101751.bcf genotype.EGAN00001101752.bcf genotype.EGAN00001101753.bcf genotype.EGAN00001101754.bcf genotype.EGAN00001101755.bcf genotype.EGAN00001101756.bcf genotype.EGAN00001101758.bcf genotype.EGAN00001101759.bcf genotype.EGAN00001101761.bcf genotype.EGAN00001101767.bcf genotype.EGAN00001101768.bcf genotype.EGAN00001101769.bcf genotype.EGAN00001101771.bcf genotype.EGAN00001101772.bcf genotype.EGAN00001101774.bcf genotype.EGAN00001101776.bcf genotype.EGAN00001101780.bcf genotype.EGAN00001101781.bcf genotype.EGAN00001101782.bcf genotype.EGAN00001101783.bcf genotype.EGAN00001101784.bcf genotype.EGAN00001101786.bcf genotype.EGAN00001101787.bcf genotype.EGAN00001101788.bcf genotype.EGAN00001101791.bcf genotype.EGAN00001101792.bcf genotype.EGAN00001101793.bcf genotype.EGAN00001101794.bcf genotype.EGAN00001101796.bcf genotype.EGAN00001101797.bcf genotype.EGAN00001101798.bcf genotype.EGAN00001101799.bcf genotype.EGAN00001101801.bcf genotype.EGAN00001101802.bcf genotype.EGAN00001101803.bcf genotype.EGAN00001101804.bcf genotype.EGAN00001101807.bcf genotype.EGAN00001101808.bcf genotype.EGAN00001101809.bcf genotype.EGAN00001101813.bcf genotype.EGAN00001101814.bcf genotype.EGAN00001101816.bcf genotype.EGAN00001101819.bcf genotype.EGAN00001101820.bcf genotype.EGAN00001101823.bcf genotype.EGAN00001101824.bcf genotype.EGAN00001101825.bcf genotype.EGAN00001101827.bcf genotype.EGAN00001101829.bcf genotype.EGAN00001101830.bcf genotype.EGAN00001101831.bcf genotype.EGAN00001101835.bcf genotype.EGAN00001101839.bcf genotype.EGAN00001101840.bcf genotype.EGAN00001101841.bcf genotype.EGYPTREF.bcf genotype.LU18.bcf genotype.LU19.bcf genotype.LU22.bcf genotype.LU23.bcf genotype.LU2.bcf genotype.LU9.bcf genotype.PD114.bcf genotype.PD115.bcf genotype.PD82.bcf
cd ..
bcftools index merged.bcf

delly filter -f germline -o germline.bcf merged.bcf
