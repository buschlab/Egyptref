#!/usr/bin/perl


use lib('/home/munz/workspace/perl_modules');


use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use File::Basename;
use POSIX;
use Cwd 'abs_path';
use Path::Class qw(dir);
use Math::Round qw(nearest);
use File::Reader;


# Directories
my $root_dir = "/data/lied_egypt_genome"; die "Directory '$root_dir' doesn't exist\n" if(!-e $root_dir);
my $programs_dir = $root_dir."/programs"; die "Directory '$programs_dir' doesn't exist\n" if(!-e $programs_dir);
my $reference_dir = $root_dir."/reference/hg38"; die "Directory '$reference_dir' doesn't exist\n" if(!-e $reference_dir);


# Programs
my $bin_java = "java";
my $bin_bwa = "bwa-0.7.17/bwa"; die "File '$bin_bwa' doesn't exist\n" if (!-e $programs_dir."/".$bin_bwa);
my $bin_tabix = "tabix-0.2.6/tabix"; die "File '$bin_tabix' doesn't exist\n" if (!-e $programs_dir."/".$bin_tabix);
my $bin_bgzip = "tabix-0.2.6/bgzip"; die "File '$bin_bgzip' doesn't exist\n" if (!-e $programs_dir."/".$bin_bgzip);
my $bin_intersect_bed = "bedtools2/bin/intersectBed"; die "File '$bin_intersect_bed' doesn't exist\n" if (!-e $programs_dir."/".$bin_intersect_bed);
my $bin_samtools = "samtools-1.3.1/samtools"; die "File '$bin_samtools' doesn't exist\n" if (!-e $programs_dir."/".$bin_samtools);
my $bin_fastqc = "FastQC/fastqc"; die "File '$bin_fastqc' doesn't exist\n" if (!-e $programs_dir."/".$bin_fastqc);
my $bin_verify_bam_id = "verifyBamID.20120620"; die "File '$bin_verify_bam_id' doesn't exist\n" if (!-e $programs_dir."/".$bin_verify_bam_id);
my $bin_vcftools = "vcftools/bin/vcftools"; die "File '$bin_vcftools' doesn't exist\n" if (!-e $programs_dir."/".$bin_vcftools);
my $pl_filterVCFbyAD = "filterVCFbyAD.pl"; die "File '$pl_filterVCFbyAD' doesn't exist\n" if (!-e $programs_dir."/".$pl_filterVCFbyAD);
my $jar_gatk_3 = "GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"; die "File '$jar_gatk_3' doesn't exist\n" if (!-e $programs_dir."/".$jar_gatk_3);
my $jar_picard = "picard.jar"; die "File '$jar_picard' doesn't exist\n" if (!-e $programs_dir."/".$jar_picard);
my $jar_trimmomatic = "Trimmomatic-0.38/trimmomatic-0.38.jar"; die "File '$jar_trimmomatic' doesn't exist\n" if (!-e $programs_dir."/".$jar_trimmomatic);
my $pl_parallelize_cmds = "parallelize_cmds.pl"; die if (!-e $programs_dir."/".$pl_parallelize_cmds);


# Reference
my $reference = "Homo_sapiens_assembly38.fasta"; die "File '$reference' doesn't exist\n" if (!-e $reference_dir."/".$reference);
my $chromosomes = "Homo_sapiens_assembly38.fasta.chr_order"; die "File '$chromosomes' doesn't exist\n" if (!-e $reference_dir."/".$chromosomes);
my $snv_dbsnp146 = "dbsnp_146.hg38.vcf.gz"; die "File '$snv_dbsnp146' doesn't exist\n" if (!-e $reference_dir."/".$snv_dbsnp146);
my $snv_hapmap = "hapmap_3.3.hg38.vcf.gz"; die "File '$snv_hapmap' doesn't exist\n" if (!-e $reference_dir."/".$snv_hapmap);
my $snv_1000GP1_hc = "1000G_phase1.snps.high_confidence.hg38.vcf.gz"; die "File '$snv_1000GP1_hc' doesn't exist\n" if (!-e $reference_dir."/".$snv_1000GP1_hc);
my $snv_1000G_omni = "1000G_omni2.5.hg38.vcf.gz"; die "File '$snv_1000G_omni' doesn't exist\n" if (!-e $reference_dir."/".$snv_1000G_omni);
my $indels_mills1000G = "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"; die "File '$indels_mills1000G' doesn't exist\n" if (!-e $reference_dir."/".$indels_mills1000G);
my $trueseq3_pe_adapters = "Trimmomatic-0.38/adapters/TruSeq4-PE.fa"; die "File '$trueseq3_pe_adapters' doesn't exist\n" if (!-e $programs_dir."/".$trueseq3_pe_adapters);
my $interval_list = "wgs_calling_regions.hg38.interval_list"; die "File '$interval_list' doesn't exist\n" if (!-e $reference_dir."/".$interval_list);
my $refseq_genes = "geneTrack.refSeq.gz"; die "File '$refseq_genes' doesn't exist\n" if (!-e $reference_dir."/".$refseq_genes);


if (!defined $ARGV[0] || !defined $ARGV[1] || !defined $ARGV[2] || !defined $ARGV[3] || !defined $ARGV[4] || !defined $ARGV[5]){

	print STDERR "Arguments required: sample_file output_root_dir mode start stop skip_base_recal\n";
	print STDERR "Level for start/stop parameters:\n1: Read group\n2: Sample\n3: Batch\n4: Study\n";
	print STDERR "Modes:\nrun: start analysis from beginning\ncontinue: continue analysis\nclean: create list of temporary files\n";
	print STDERR "Set skip_base_recal_stats=1, if calculation of base recalibration statistics should be skipped, otherwise skip_base_recal_stats=0\n";

	exit(0);
}


my ($sample_file, $output_root_dir, $mode, $start, $stop, $skip_base_recal_stats) = @ARGV;


die "File '$sample_file' doesn't exist\n" if (!-e $sample_file);
die "Unknown mode '$mode'\n" if ($mode ne "run" && $mode ne "continue" && $mode ne "clean" && $mode ne "error");
if(!looks_like_number($start) || !looks_like_number($stop)){
	$start = 1;
	$stop = 4;
}
die "Unknown value for skip_base_recal: '$skip_base_recal_stats'\n" if (!looks_like_number($skip_base_recal_stats) or ($skip_base_recal_stats != 1 && $skip_base_recal_stats != 0));


my %chromosomes;
my @chromosomes;
open(IN, "<$reference_dir/$chromosomes") or die "Can't open file '$reference_dir/$chromosomes': $!\n";
while(<IN>){
	chomp($_);
	die if(exists($chromosomes{$_}));
	$chromosomes{$_} = 1;
	push(@chromosomes, $_);
}


print "Sample File: $sample_file\n";
print "Output Root Dir: $output_root_dir\n";
print "Mode: $mode\n";
print "Start: $start\n";
print "Stop: $stop\n";
print "Skip base recalibration statistics: $skip_base_recal_stats\n";
print "Chromosomes: ".join(" ", @chromosomes)."\n";


# Create output root dir
$output_root_dir = abs_path($output_root_dir);
system("mkdir -p $output_root_dir");
die "Error: Can't create '$output_root_dir'\n" if (!-e $output_root_dir);


# Tmp files
my @tmp_files;


# Translate mode
my $continue = 0;
if ($mode eq "continue"){
	$continue = 1;
}


# -----------------------------------------------
# Read sample data
# -----------------------------------------------
print "...read '$sample_file'\n";
my %sample_id2rg2normrg_fastq1_fastq2;
my $n_samples = 0;
my $n_read_groups = 0;
my $reader = File::Reader -> new($sample_file, {'has_header' => 1, 'skip' => '#'});
my ($header, $header_inv, $header_mult) = $reader -> get_header();
my $sample_colname = "sample_id"; die if(!exists(${$header}{$sample_colname}));
my $lib_colname = "lib"; die if(!exists(${$header}{$lib_colname}));
my $flowcell_colname = "flowcell"; die if(!exists(${$header}{$flowcell_colname}));
my $lane_colname = "lane"; die if(!exists(${$header}{$lane_colname}));
my $barcode_colname = "barcode"; die if(!exists(${$header}{$barcode_colname}));
my $platform_colname = "platform"; die if(!exists(${$header}{$platform_colname}));
my $fastq1_colname = "fastq1"; die if(!exists(${$header}{$fastq1_colname}));
my $fastq2_colname = "fastq2"; die if(!exists(${$header}{$fastq2_colname}));
while($reader -> has_next()){

		my @row = $reader -> next();
		next if($row[0] =~ /^#/);

		my ($sample_id, $flowcell, $lane, $barcode, $lib, $platform, $fastq1, $fastq2) = ($row[${$header}{$sample_colname}], $row[${$header}{$flowcell_colname}], $row[${$header}{$lane_colname}], $row[${$header}{$barcode_colname}], $row[${$header}{$lib_colname}], $row[${$header}{$platform_colname}], $row[${$header}{$fastq1_colname}], $row[${$header}{$fastq2_colname}]);

		die if (!defined $sample_id || $sample_id eq "");
		die if (!defined $flowcell || $flowcell eq "");
		die if (!defined $lane || $lane eq "");
		die if (!defined $barcode || $barcode eq "");
		die if (!defined $lib || $lib eq "");
		die "File '$fastq1' doesn't exist\n" if (!-e $fastq1);
		die "File '$fastq2' doesn't exist\n" if (!-e $fastq2);
		die "Files are the same: '$fastq1' vs '$fastq2'\n" if($fastq1 eq $fastq2);


		# LB:<sample_id>_<lib>
		# PL: Illumina
		# ID: <FLOWCELL>_<LANE>
		# PU: <FLOWCELL>_<LANE>_<SAMPLE_BARCODE>
		my $rg_id = "${flowcell}_$lane";
		my $pu = "${flowcell}_${lane}_$barcode";
		my $l = "${sample_id}_$lib";
		my $read_group = join("\\t", ("ID:$rg_id", "SM:$sample_id", "LB:$l", "PU:$pu", "PL:$platform"));
		my $norm_rg = join("_", ($lib, $flowcell, $lane, $barcode));
		$norm_rg =~ s/[^A-Za-z0-9]/_/g; # replace all non-alphanumericals with "_"
		$norm_rg =~ s/_+/_/g;
		print "Read Group: ".join(" ", ("ID:$rg_id", "SM:$sample_id", "LB:$lib", "PU:$pu", "PL:$platform"))." (Normalized: $norm_rg)\n";

		die if (exists($sample_id2rg2normrg_fastq1_fastq2{$sample_id}) && exists($sample_id2rg2normrg_fastq1_fastq2{$sample_id}{$read_group}));

		$n_samples++ if (!exists($sample_id2rg2normrg_fastq1_fastq2{$sample_id}));
		$n_read_groups++;

		$sample_id2rg2normrg_fastq1_fastq2{$sample_id}{$read_group} = [$norm_rg, $fastq1, $fastq2];

}
print "\t#samples: $n_samples\n";
print "\t#read groups: $n_read_groups\n";


# -----------------------------------------------
# Processing per sample read group
# -----------------------------------------------
print "...create sample read group processing scripts\n";
my %sample_id2job_ids;
my %sample_id2mapping_files;
my %sample_id2continue_rg;
my $level = 1;
foreach my $sample_id (sort (keys(%sample_id2rg2normrg_fastq1_fastq2))){
	foreach my $read_group (sort (keys(%{$sample_id2rg2normrg_fastq1_fastq2{$sample_id}}))){

		my $continue_rg = $continue;

		my ($norm_rg, $fastq_r1, $fastq_r2) = @{$sample_id2rg2normrg_fastq1_fastq2{$sample_id}{$read_group}};
		my ($output_dir) = _create_output_dir([$sample_id, $norm_rg], $output_root_dir);

		my $job_id;


		# FastQC
		if($start <= $level && $stop >= $level){

			my $fastqc_dir = "fastqc_01_raw";

			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my $cores = 22;
				my $mem = "20GB";

				#TODO: Check scratch size before & after

				my $time = 700;

				my @cmds;
				push(@cmds, _cp_cmd($fastq_r1, $fastq_r2, '$SCRATCH'));
				push(@cmds, _fastqc_cmd(_path($programs_dir, $bin_fastqc), $cores, _path('$SCRATCH', $fastqc_dir, basename($fastq_r1), basename($fastq_r2))));
				push(@cmds, _cp_cmd(_path('$SCRATCH', $fastqc_dir), $output_dir));

				my $sh_file = "$output_dir/$sample_id.$norm_rg.$fastqc_dir.sh";
				my $log_file = "$output_dir/$sample_id.$norm_rg.$fastqc_dir.log";
				(undef, undef) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $fastqc_dir)], $continue_rg, $cores, $mem, $time, [], $mode);
			}
		}


		my $trimmomatic_forward_paired_file = "$sample_id.$norm_rg.R1.trimmomatic.fq.gz";
		my $trimmomatic_reverse_paired_file = "$sample_id.$norm_rg.R2.trimmomatic.fq.gz";
		if($start <= $level && $stop >= $level){

			my $trimmomatic_forward_unpaired_file = "$sample_id.$norm_rg.R1.trimmomatic.unpaired.fq.gz";
			my $trimmomatic_reverse_unpaired_file = "$sample_id.$norm_rg.R2.trimmomatic.unpaired.fq.gz";
			my $trimmomatic_log_file = "$sample_id.$norm_rg.trimmomatic.log";

			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my $cores = 10;
				my $mem = "10GB";  my $mem_java = "8g";
				my $time = 1000;

				my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
				push(@cmds, _cp_cmd($fastq_r1, $fastq_r2, '$SCRATCH'));
				push(@cmds, _trimmomatic_cmd($bin_java, _path("\$SCRATCH/programs", $jar_trimmomatic, $trueseq3_pe_adapters), $cores, $mem_java, _path('$SCRATCH', basename($fastq_r1), basename($fastq_r2), $trimmomatic_forward_paired_file, $trimmomatic_reverse_paired_file, $trimmomatic_forward_unpaired_file, $trimmomatic_reverse_unpaired_file, $trimmomatic_log_file)));
				push(@cmds, _cp_cmd(_path('$SCRATCH', $trimmomatic_forward_paired_file, $trimmomatic_reverse_paired_file), $output_dir));

				my $sh_file = "$output_dir/$sample_id.$norm_rg.trimmomatic.sh";
				my $log_file = "$output_dir/$sample_id.$norm_rg.trimmomatic.log";
				($job_id, $continue_rg) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $trimmomatic_forward_paired_file, $trimmomatic_reverse_paired_file)], $continue, $cores, $mem, $time, [], $mode);
			}


			push(@tmp_files, _path($output_dir, $trimmomatic_forward_paired_file, $trimmomatic_reverse_paired_file));
		}


		# FastQC
		if($start <= $level && $stop >= $level){

			my $fastqc_dir = "fastqc_02_trimmomatic";

			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my $cores = 22;
				my $mem = "20GB";
				my $time = 700;

				my @cmds;
				push(@cmds, _cp_cmd($fastq_r1, $fastq_r2, '$SCRATCH'));
				push(@cmds, _fastqc_cmd(_path($programs_dir, $bin_fastqc), $cores, _path('$SCRATCH', $fastqc_dir, basename($fastq_r1), basename($fastq_r2))));
				push(@cmds, _cp_cmd(_path('$SCRATCH', $fastqc_dir), $output_dir));

				my $sh_file = "$output_dir/$sample_id.$norm_rg.$fastqc_dir.sh";
				my $log_file = "$output_dir/$sample_id.$norm_rg.$fastqc_dir.log";
				(undef, undef) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $fastqc_dir)], $continue_rg, $cores, $mem, $time, [$job_id], $mode);
			}
		}


		# Mapping to reference using BWA mem
		my $mapping_file = "$sample_id.$norm_rg.mapping.bam";
		if($start <= $level && $stop >= $level){


			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my $cores = 14;
				my $mem = "70GB"; my $mem_samtools = "4g";
				my $time = 2000;

				my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
				push(@cmds, _cp_cmd(_path($output_dir, $trimmomatic_forward_paired_file, $trimmomatic_reverse_paired_file), '$SCRATCH'));
				push(@cmds, _mapping_cmd($read_group, _path('$SCRATCH/programs', $bin_bwa, $bin_samtools), _path('$SCRATCH/reference', $reference), _path('$SCRATCH', $trimmomatic_forward_paired_file, $trimmomatic_reverse_paired_file, $mapping_file), $cores, $mem_samtools));
				push(@cmds, _bam_index_cmd(_path('$SCRATCH/programs', $bin_samtools), _path('$SCRATCH', $mapping_file)));
				push(@cmds, _cp_cmd(_path('$SCRATCH', $mapping_file, $mapping_file.".bai"), $output_dir));

				my $sh_file = "$output_dir/$sample_id.$norm_rg.mapping.sh";
				my $log_file = "$output_dir/$sample_id.$norm_rg.mapping.log";
				($job_id, $continue_rg) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $mapping_file, $mapping_file.".bai")], $continue_rg, $cores, $mem, $time, [$job_id], $mode);


				push(@{$sample_id2job_ids{$sample_id}}, $job_id) if(looks_like_number($job_id));
			}


			push(@{$sample_id2mapping_files{$sample_id}}, "$output_dir/$mapping_file");
			push(@tmp_files, _path($output_dir, $mapping_file));
		}
		else {
			push(@{$sample_id2mapping_files{$sample_id}}, "$output_dir/$mapping_file");
		}


		# Continue
		if(!exists($sample_id2continue_rg{$sample_id})){
			$sample_id2continue_rg{$sample_id} = $continue_rg;
		}
		elsif($sample_id2continue_rg{$sample_id} != 0){
			$sample_id2continue_rg{$sample_id} = $continue_rg;
		}
	}
}

if ($stop <= 1){
	_create_tmp_file_list(\@tmp_files, $output_root_dir."/tmp_files.txt") if($mode eq "clean");
	exit;
}


# -----------------------------------------------
# Processing per sample
# -----------------------------------------------
print "...create sample scripts\n";
my %gvcf_file2job_id;
my %gvcf_file2continue;
my %sample_id2bam_file;
my %bam_file2job_id;
my %bam_file2continue;
$level = 2;
foreach my $sample_id (sort (keys(%sample_id2rg2normrg_fastq1_fastq2))){

	my @rg_job_ids = @{$sample_id2job_ids{$sample_id}} if (exists($sample_id2job_ids{$sample_id}));
	my $continue_sample = $sample_id2continue_rg{$sample_id};

	my ($output_dir) = _create_output_dir([$sample_id], $output_root_dir);

	my $job_id;


	# Merge BAMs
	my $merged_mapping_file = "$sample_id.merged.bam";
	my $merged_mapping_file_dir;
	if($start <= $level && $stop >= $level){
		my @mapping_files = @{$sample_id2mapping_files{$sample_id}};

		if (@mapping_files > 1){

			$merged_mapping_file_dir = $output_dir;


			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my $cores = 1;
				my $mem = "10GB";
				my $time = 1200;

				my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
				push(@cmds, _cp_cmd(@mapping_files, '$SCRATCH'));
				push(@cmds, _merge_bam_files_cmd(_path('$SCRATCH/programs', $bin_samtools), _path('$SCRATCH', $merged_mapping_file, _basename_list(@mapping_files))));
				push(@cmds, _bam_index_cmd(_path('$SCRATCH/programs', $bin_samtools), _path('$SCRATCH', $merged_mapping_file)));
				push(@cmds, _cp_cmd(_path('$SCRATCH', $merged_mapping_file, $merged_mapping_file.".bai"), $output_dir));

				my $sh_file = "$output_dir/$sample_id.merged.sh";
				my $log_file = "$output_dir/$sample_id.merged.log";
				($job_id, $continue_sample) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $merged_mapping_file, $merged_mapping_file.".bai")], $continue_sample, $cores, $mem, $time, [@rg_job_ids], $mode);
			}


			push(@tmp_files, _path($output_dir, $merged_mapping_file));
		}
		elsif(@mapping_files == 1){

			$merged_mapping_file = basename($mapping_files[0]);
			$merged_mapping_file_dir = dirname($mapping_files[0]);
		}
		else {
			die;
		}
	}


#	# Blast unmapped reads
#	$step_id_sample++;
#	if($start <= $level && $stop >= $level){
#
#		my ($output_dir, '$SCRATCH', '$SCRATCH/reference', '$SCRATCH/programs') = _create_project_dirs([$sample_id], $step_id_sample, $output_root_dir);
#
#		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){
#
#			my $cores = 1;
#			my $mem = "1GB";
#			my $time = 400;
#
#
#			my $sample_size = 1000;
#			my $unmapped_reads_file = "unmapped_reads.n$sample_size.fa";
#			my $blastn_file = "unmapped_reads.n$sample_size.blast.txt";
#
#
#			my @cmds;
#			push(@cmds, "sleep ".int(rand($time-30))."m");
#			push(@cmds, _time(_blast_unmapped_reads(_path($programs_dir, $bin_samtools), "/home/munz/workspace/blast/web_blast.pl",  _path($merged_mapping_file_dir, $merged_mapping_file),_path($output_dir, $unmapped_reads_file, $blastn_file), $sample_size)));
#			push(@cmds, "echo '\@COMPLETE'");
#
#
#			my $sh_file = "$output_dir/$sample_id.blast.sh";
#			my $log_file = "$output_dir/$sample_id.blast.log";
#			(undef, undef) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $unmapped_reads_file, $blastn_file)], $continue_sample, $cores, $mem, $time, [], $mode);
#
#		}
#	}


	# Mark duplicates
	my $marked_duplicates_file = "$sample_id.merged.mark_dups.bam";
	if($start <= $level && $stop >= $level){


		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $duplication_metrics_file = "$sample_id.marked_dup_metrics.txt";

			my $cores = 4;
			my $mem = "80GB"; my $mem_java = "75g";
			my $time = 1000;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs");
			push(@cmds, _cp_cmd(_path($merged_mapping_file_dir, $merged_mapping_file, $merged_mapping_file.".bai"), '$SCRATCH'));
			push(@cmds, _mark_duplicates($bin_java, $mem_java, $cores, _path('$SCRATCH/programs', $jar_picard), _path('$SCRATCH', $merged_mapping_file, $marked_duplicates_file, $duplication_metrics_file), '$SCRATCH'));
			push(@cmds, _bam_index_cmd(_path('$SCRATCH/programs', $bin_samtools), _path('$SCRATCH', $marked_duplicates_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $duplication_metrics_file, $marked_duplicates_file, $marked_duplicates_file.".bai"), $output_dir));

			my $sh_file = "$output_dir/$sample_id.mark_dups.sh";
			my $log_file = "$output_dir/$sample_id.mark_dups.log";
			($job_id, $continue_sample) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $marked_duplicates_file, $marked_duplicates_file.".bai", $duplication_metrics_file)], $continue_sample, $cores, $mem, $time, [$job_id, @rg_job_ids], $mode);

		}
	}


	# Base recalibration
	my $recal_report_file = "$sample_id.base_recal.table";
	if($start <= $level && $stop >= $level){

		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 12;
			my $mem = "80GB"; my $mem_java = "75g";
			my $time = 2000;

			my $known_sites = {_path('$SCRATCH/reference', $snv_dbsnp146) => 1, _path('$SCRATCH/reference',  $indels_mills1000G) => 1};

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $marked_duplicates_file, $marked_duplicates_file.".bai"), '$SCRATCH'));
			push(@cmds, _base_recal_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), $cores, $mem_java, _path('$SCRATCH/reference', $reference), $known_sites, _path('$SCRATCH', $marked_duplicates_file, $recal_report_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $recal_report_file), $output_dir));

			my $sh_file = "$output_dir/$sample_id.base_recal.sh";
			my $log_file = "$output_dir/$sample_id.base_recal.log";
			($job_id, $continue_sample) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $recal_report_file)], $continue_sample, $cores, $mem, $time, [$job_id, @rg_job_ids], $mode);
		}
	}


	# Apply base recalibration
	my $recal_bam_file = "$sample_id.merged.mark_dups.base_recal.bam";
	$sample_id2bam_file{$sample_id} = "$output_dir/$recal_bam_file";
	$bam_file2continue{"$output_dir/$recal_bam_file"} = "NA";
	$bam_file2job_id{"$output_dir/$recal_bam_file"} = "NA";
	if($start <= $level && $stop >= $level){

		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 12;
			my $mem = "50GB"; my $mem_java = "45g";
			my $time = 1800;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $marked_duplicates_file, $marked_duplicates_file.".bai", $recal_report_file), '$SCRATCH'));
			push(@cmds, _apply_base_recal_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), $cores, $mem_java, _path('$SCRATCH/reference', $reference), _path('$SCRATCH', $marked_duplicates_file, $recal_report_file, $recal_bam_file)));
			push(@cmds, _bam_index_cmd(_path('$SCRATCH/programs', $bin_samtools), _path('$SCRATCH', $recal_bam_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $recal_bam_file, $recal_bam_file.".bai"), $output_dir));

			my $sh_file = "$output_dir/$sample_id.apply_base_recal.sh";
			my $log_file = "$output_dir/$sample_id.apply_base_recal.log";
			($job_id, $continue_sample) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $recal_bam_file, $recal_bam_file.".bai")], $continue_sample, $cores, $mem, $time, [$job_id], $mode);

			$bam_file2continue{"$output_dir/$recal_bam_file"} = $continue_sample;
			$bam_file2job_id{"$output_dir/$recal_bam_file"} = $job_id;
		}
	}


	# Base recalibration statistics
	if(!$skip_base_recal_stats && $start <= $level && $stop >= $level){

		my $post_recal_report_file = "$sample_id.base_recal.post.table";
		my $recal_plots = "$sample_id.base_recal.plots.pdf";

		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 12;
			my $mem = "80GB"; my $mem_java = "75g";
			my $time = 8000;

			my $known_sites = {_path('$SCRATCH/reference', $snv_dbsnp146) => 1, _path('$SCRATCH/reference',  $indels_mills1000G) => 1};

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $marked_duplicates_file, $marked_duplicates_file.".bai", $recal_report_file), '$SCRATCH'));
			push(@cmds, _base_recal_stats_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), $cores, $mem_java, _path('$SCRATCH/reference', $reference), $known_sites, _path('$SCRATCH', $marked_duplicates_file, $recal_report_file, $post_recal_report_file, $recal_plots)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $post_recal_report_file, $recal_plots), $output_dir));

			my $sh_file = "$output_dir/$sample_id.base_recal_stats.sh";
			my $log_file = "$output_dir/$sample_id.base_recal_stats.log";
			(undef, undef) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $post_recal_report_file, $recal_plots)], $continue_sample, $cores, $mem, $time, [$job_id], $mode);
		}
	}


	# Picard metrics
	if($start <= $level && $stop >= $level){

		my $picard_prefix = "picard/$sample_id";

		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 9;
			my $mem = "95GB"; my $mem_java = "10g";
			my $time = 1200;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference", "mkdir -p \$SCRATCH/".dirname($picard_prefix));
			push(@cmds, _cp_cmd(_path($output_dir, $recal_bam_file, $recal_bam_file.".bai"), '$SCRATCH'));
			push(@cmds, _picard_metrics_cmd($bin_java, _path('$SCRATCH/programs', $jar_picard), $mem_java, _path('$SCRATCH/reference', $reference), _path('$SCRATCH', $recal_bam_file, $picard_prefix)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', dirname($picard_prefix)), $output_dir));

			my $sh_file = "$output_dir/$sample_id.picard.sh";
			my $log_file = "$output_dir/$sample_id.picard.log";
			(undef, undef) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, dirname($picard_prefix))], $continue_sample, $cores, $mem, $time, [$job_id], $mode);
		}
	}


	# Per sample HaplotypeCaller
	my $gvcf_file = "$sample_id.g.vcf.gz";
	if($start <= $level && $stop >= $level){

		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 8;
			my $mem = "30GB"; my $mem_java = "28g";
			my $time = 4000;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $recal_bam_file, $recal_bam_file.".bai"), '$SCRATCH'));
			push(@cmds, _hc_per_sample_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), _path('$SCRATCH/reference', $reference, $interval_list), $cores, $mem_java, _path('$SCRATCH', $recal_bam_file, $gvcf_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $gvcf_file, $gvcf_file.".tbi"), $output_dir));

			my $sh_file = "$output_dir/$sample_id.hc.sh";
			my $log_file = "$output_dir/$sample_id.hc.log";
			($job_id, $continue_sample) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $gvcf_file, $gvcf_file.".tbi")], $continue_sample, $cores, $mem, $time, [$job_id], $mode);

			$gvcf_file2job_id{$output_dir."/".$gvcf_file} = $job_id;
			$gvcf_file2continue{$output_dir."/".$gvcf_file} = $continue_sample;
		}
		else {
			$gvcf_file2job_id{$output_dir."/".$gvcf_file} = "NA";
			$gvcf_file2continue{$output_dir."/".$gvcf_file} = $continue_sample;
		}
	}
	else {
		$gvcf_file2job_id{$output_dir."/".$gvcf_file} = "NA";
		$gvcf_file2continue{$output_dir."/".$gvcf_file} = $continue_sample;
	}
}

if ($stop <= 2){
	_create_tmp_file_list(\@tmp_files, $output_root_dir."/tmp_files.txt") if($mode eq "clean");
	exit;
}


$level = 3;
{

	# -----------------------------------------------
	# Coverage depth
	# -----------------------------------------------
	if($start <= $level && $stop >= $level){

		my $output_dir = $output_root_dir."/cov_depth";
		system("mkdir -p $output_dir"); die if (!-e $output_dir);


		my $covdepth_list_file = $output_dir."/bam_files.list";
		open(OUT, ">".$covdepth_list_file) or die "Cannot open file '$covdepth_list_file': $!\n";
		print OUT join("\n", sort keys(%bam_file2continue))."\n";
		close OUT;

#TODO: Merge bams to scratch, create index and run per chromosomes
#		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){
#
#
#			my @covdepth_depend; push(@covdepth_depend, $bam_file2job_id{$_}) for keys(%bam_file2job_id);
#			my $covdepth_continue = 1;
#			foreach my $bam_file (keys(%bam_file2continue)){
#				if(looks_like_number($bam_file2continue{$bam_file}) && $bam_file2continue{$bam_file} == 0){
#					$covdepth_continue = 0;
#				}
#			}
#
#			my $cores = 1;
#			my $mem = "8GB"; my $mem_java = "8g";
#			my $time = 7200;
#
#			my @cmds;
#			push(@cmds, _cov_depth_cmd($bin_java, _path($programs_dir, $jar_gatk_3), $cores, $mem_java, _path($reference_dir, $reference, $refseq_genes), $covdepth_list_file, $output_dir."/cov_depth"));
#
#			my $sh_file = "$output_dir/cov_depth.sh";
#			my $log_file = "$output_dir/cov_depth.log";
#			(undef, undef) = _run(\@cmds, $sh_file, $log_file, [], $covdepth_continue, $cores, $mem, $time, [@covdepth_depend], $mode);
#		}
	}
}


my %merged_gvcf_file2job_id;
my $continue_study = 1;
{
	# -----------------------------------------------
	# Create batch lists
	# -----------------------------------------------
	my $min_gvcfs_per_batch = 50;
	my $max_gvcfs_per_batch = 100;

	print "...create batch lists\n";
	print "\tmax number of samples per batch: $max_gvcfs_per_batch\n";
	print "\tmin number of samples per batch: $min_gvcfs_per_batch\n";


	my @batches = _get_best_split([sort keys(%gvcf_file2job_id)], $min_gvcfs_per_batch, $max_gvcfs_per_batch);
	my %batch_id2continue;
	my %batch_id2gvcf_file2job_id;
	for (my $i=1; $i<=@batches; $i++){

		my $batch_id = "batch_$i";

		my $batch_dir = $output_root_dir."/$batch_id";
		system("mkdir -p $batch_dir"); die if (!-e $batch_dir);

		my $batch_list_file = $batch_dir."/${batch_id}.gvcfs.txt";
		open(OUT, ">".$batch_list_file) or die "Cannot open file '$batch_list_file': $!\n";
		print OUT join("\n", @{$batches[$i-1]})."\n";
		close OUT;


		foreach my $gvcf (@{$batches[$i-1]}){
			$batch_id2gvcf_file2job_id{$batch_id}{$gvcf} = $gvcf_file2job_id{$gvcf};


			# Continue
			if(!exists($batch_id2continue{$batch_id}) || $batch_id2continue{$batch_id} != 0){
				$batch_id2continue{$batch_id} = $gvcf_file2continue{$gvcf};
			}
		}
	}
	my $n_batches = @batches;
	print "\t$n_samples samples were distributed over $n_batches batches\n";


	# -----------------------------------------------
	# Merge to batches
	# -----------------------------------------------
	print "...create batch scripts\n";
	foreach my $batch_id (sort keys(%batch_id2gvcf_file2job_id)){

		my $continue_batch = $batch_id2continue{$batch_id};
		my ($output_dir) = _create_output_dir([$batch_id], $output_root_dir);
		my $job_id;


		# Merge gvcf files (max. 100 per batch)
		my $batch_merge_file_prefix = $batch_id;
		if($start <= $level && $stop >= $level){

			my @sample_job_ids;
			my @gvcf_files;
			foreach my $gvcf_file (sort keys(%{$batch_id2gvcf_file2job_id{$batch_id}})){
				push(@sample_job_ids, $batch_id2gvcf_file2job_id{$batch_id}{$gvcf_file});
				push(@gvcf_files, $gvcf_file);
			}


			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my @gvcf_files_basename;
				push(@gvcf_files_basename, basename($_)) for @gvcf_files;

				my $cores = 40;
				my $mem = "100GB"; my $mem_java = "4g";
				my $time = 7200;

				my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
				push(@cmds, _cp_cmd($_, $_.".tbi", '$SCRATCH')) for @gvcf_files;
				push(@cmds, "cd $output_dir"); # Otherwise the -L argument in SelectVariants can be recognized as sample file and not as chr
				push(@cmds, _parallel_merge_gvcfs_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3, $pl_parallelize_cmds),_path('$SCRATCH/reference', $reference), \@chromosomes, $cores, $mem_java, [_path('$SCRATCH', @gvcf_files_basename)], _path('$SCRATCH', $batch_merge_file_prefix)));
				push(@cmds, _cp_cmd(_path('$SCRATCH', $batch_merge_file_prefix.'.g.vcf.gz', $batch_merge_file_prefix.'.g.vcf.gz.tbi'), $output_dir));

				my $sh_file = "$output_dir/$batch_id.merge.sh";
				my $log_file = "$output_dir/$batch_id.merge.log";
				($job_id, $continue_batch) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $batch_merge_file_prefix.".g.vcf.gz", $batch_merge_file_prefix.".g.vcf.gz.tbi")], $continue_batch, $cores, $mem, $time, [@sample_job_ids], $mode);
			}
		}


		if($continue_batch == 0){
			$continue_study = $continue_batch;
		}


		$merged_gvcf_file2job_id{$_} = $job_id for _path($output_dir, $batch_merge_file_prefix.".g.vcf.gz");
	}
}


# -----------------------------------------------
# Processing per study
# -----------------------------------------------
$level = 4;
{

	my ($output_dir) = _create_output_dir([], $output_root_dir);
	my $job_id;


	# Multi sample variant calling
	my $hc_file = "all.hc.vcf.gz";
	if($start <= $level && $stop >= $level){


		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 4;
			my $mem = "32GB"; my $mem_java = "32g";
			my $time = 3000;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _multi_sample_hc_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), _path('$SCRATCH/reference', $reference, $snv_dbsnp146), $cores, $mem_java, [keys(%merged_gvcf_file2job_id)], _path('$SCRATCH', $hc_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH',  $hc_file, $hc_file.".tbi"), $output_dir));

			my $sh_file = "$output_dir/hc.sh";
			my $log_file = "$output_dir/hc.log";

			my @job_dependencies; push(@job_dependencies, $merged_gvcf_file2job_id{$_}) for keys(%merged_gvcf_file2job_id);
			($job_id, $continue_study) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $hc_file, $hc_file.".tbi")], $continue_study, $cores, $mem, $time, \@job_dependencies, $mode);
		}
	}


	# Variant Quality Score Recalibration (SNVs)
	my $snv_file = "all.recal_snv.vcf.gz";
	if($start <= $level && $stop >= $level){

		my $recal_snv_file = "all.recal_snv.recal";
		my $tranches_snv_file = "all.recal_snv.tranches";
		my $rscript_file = "all.recal_snv.plots.R";

		my $resources = {"hapmap" => {"known" => "false", "training" => "true", "truth" => "true", "prior" => "15.0", "file" => _path('$SCRATCH/reference', $snv_hapmap)},
						 "omni" => {"known" => "false", "training" => "true", "truth" => "false", "prior" => "12.0", "file" => _path('$SCRATCH/reference', $snv_1000G_omni)},
						 "1000G" => {"known" => "false", "training" => "true", "truth" => "false", "prior" => "10.0", "file" => _path('$SCRATCH/reference', $snv_1000GP1_hc)},
						 "dbsnp" => {"known" => "true", "training" => "false", "truth" => "false", "prior" => "2.0", "file" => _path('$SCRATCH/reference', $snv_dbsnp146)},
		};

		if(keys(%{$resources}) == 0){
			$snv_file = $hc_file;
		}
		elsif($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 4;
			my $mem = "18GB"; my $mem_java = "16g";
			my $time = 2000;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $hc_file, $hc_file.".tbi"), '$SCRATCH'));
			push(@cmds, _recal_snv_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), $cores, $mem_java, _path('$SCRATCH/reference', $reference), $resources, _path('$SCRATCH', $hc_file, $recal_snv_file, $tranches_snv_file, $rscript_file, $snv_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $snv_file, $snv_file.".tbi", $recal_snv_file, $tranches_snv_file, $tranches_snv_file.".pdf", $rscript_file, $rscript_file.".pdf"), $output_dir));

			my $sh_file = "$output_dir/all.recal_snv.sh";
			my $log_file = "$output_dir/all.recal_snv.log";
			($job_id, $continue_study) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $snv_file, $snv_file.".tbi", $recal_snv_file, $tranches_snv_file, $tranches_snv_file.".pdf", $rscript_file, $rscript_file.".pdf")], $continue_study, $cores, $mem, $time, [$job_id], $mode);
		}
	}


	# Variant Quality Score Recalibration (INDELs)
	my $recal_vcf_file = "all.recal_snv.recal_indel.vcf.gz";
	if($start <= $level && $stop >= $level){

		my $recal_indel_file = "all.recal_indel.recal";
		my $tranches_indels_file = "all.recal_indel.tranches";
		my $rscript_file = "all.recal_indel.plots.R";

		my $resources = {"mills" => {"known" => "false", "training" => "true", "truth" => "true", "prior" => "12.0", "file" => _path('$SCRATCH/reference', $indels_mills1000G)},
						 "dbsnp" => {"known" => "true", "training" => "false", "truth" => "false", "prior" => "2.0", "file" => _path('$SCRATCH/reference', $snv_dbsnp146)},
		};


		if(keys(%{$resources}) == 0){
			$recal_vcf_file = $snv_file;
		}
		elsif($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 4;
			my $mem = "18GB"; my $mem_java = "46g";
			my $time = 2000;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $snv_file, $snv_file.".tbi"), '$SCRATCH'));
			push(@cmds, _recal_indel_cmd($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), $cores, $mem_java, _path('$SCRATCH/reference', $reference), $resources, _path('$SCRATCH', $snv_file, $recal_indel_file, $tranches_indels_file, $rscript_file, $recal_vcf_file)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $recal_indel_file, $recal_vcf_file, $recal_vcf_file.".tbi", $tranches_indels_file, $rscript_file, $rscript_file.".pdf"), $output_dir));

			my $sh_file = "$output_dir/all.recal_indel.sh";
			my $log_file = "$output_dir/all.recal_indel.log";
			($job_id, $continue_study) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $recal_vcf_file, $recal_vcf_file.".tbi", $recal_indel_file, $tranches_indels_file, $rscript_file, $rscript_file.".pdf")], $continue_study, $cores, $mem, $time, [$job_id], $mode);
		}
	}


	# Create clean variant sets
	my $clean_snvs = "snv.clean.vcf.gz";
	if($start <= $level && $stop >= $level){

		my $clean_vars = "vars.clean.vcf.gz";
		my $clean_indels = "indels.clean.vcf.gz";

		if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

			my $cores = 1;
			my $mem = "5GB"; my $mem_java = "4g";
			my $time = 3000;

			my @cmds = ("cp -rf $programs_dir \$SCRATCH/programs", "cp -rf $reference_dir \$SCRATCH/reference");
			push(@cmds, _cp_cmd(_path($output_dir, $recal_vcf_file, $recal_vcf_file.".tbi"), '$SCRATCH'));
			push(@cmds, _create_clean_variant_sets($bin_java, _path('$SCRATCH/programs', $jar_gatk_3), $mem_java, _path('$SCRATCH/reference', $reference), _path('$SCRATCH', $recal_vcf_file, $clean_vars, $clean_snvs, $clean_indels)));
			push(@cmds, _cp_cmd(_path('$SCRATCH', $clean_snvs, $clean_snvs.".tbi", $clean_indels, $clean_indels.".tbi", $clean_vars, $clean_vars.".tbi"), $output_dir));

			my $sh_file = "$output_dir/clean.sh";
			my $log_file = "$output_dir/clean.log";
			($job_id, $continue_study) = _run(\@cmds, $sh_file, $log_file, [_path($output_dir, $clean_snvs, $clean_snvs.".tbi", $clean_indels, $clean_indels.".tbi", $clean_vars, $clean_vars.".tbi")], $continue_study, $cores, $mem, $time, [$job_id], $mode);
		}
	}


	foreach my $sample_id (sort (keys(%sample_id2rg2normrg_fastq1_fastq2))){

		# VerifyBamID
		if($start <= $level && $stop >= $level){

			my ($sample_output_dir) = _create_output_dir([$sample_id], $output_root_dir);
			my $verify_bam_id_output_prefix = "$sample_id.verifybamid";

			if($mode eq "run" or $mode eq "continue" or $mode eq "error"){

				my $cores = 1;
				my $mem = "10GB";
				my $time = 5000;

				my @cmds;
				push(@cmds, _verify_bam_id_cmd(_path($programs_dir, $bin_verify_bam_id), _path($output_dir, $clean_snvs), $sample_id2bam_file{$sample_id}, _path($sample_output_dir, $verify_bam_id_output_prefix)));

				my $sh_file = "$sample_output_dir/$sample_id.verifybamid.sh";
				my $log_file = "$sample_output_dir/$sample_id.verifybamid.log";
				(undef, undef) = _run(\@cmds, $sh_file, $log_file, [_path($sample_output_dir, "$sample_id.verifybamid.selfSM", "$sample_id.verifybamid.depthSM", "$sample_id.verifybamid.bestSM", "$sample_id.verifybamid.bestRG", "$sample_id.verifybamid.selfRG", "$sample_id.verifybamid.depthRG")], $continue_study, $cores, $mem, $time, [$job_id], $mode);
			}
		}
	}
}


_create_tmp_file_list(\@tmp_files, $output_root_dir."/tmp_files.txt") if($mode eq "clean");


sub _create_tmp_file_list {
	my ($files, $output_file) = @_;

	# -----------------------------------------------
	# Create list with temporary files
	# -----------------------------------------------

	print "...write temporary files to '$output_file'\n";
	open(OUT, ">".$output_file) or die "Cannot open file: $!\n";
	my $sum = 0;
	foreach my $f (@{$files}){

		my ($s) = $f =~ /^(.*?)\**$/;

		if(!-e $s){
			print STDERR "File '$f' doesn't exist\n";
		}
		else {
			my $size_in_bytes = -s $s;
			$sum += $size_in_bytes;
		}

		print OUT $f."\n";
	}
	close OUT;

	my $size_in_tb = nearest(.1, $sum/(1024**4));

	print "Overall size of temporary files: ".$size_in_tb."TB\n";
	print "Type 'xargs -I\@ sh -c 'rm -rf \@' <$output_file' to delete temporary files\n";
}


# -----------------------------------------------
sub _fastqc_cmd {
# -----------------------------------------------
	my ($bin_fastqc, $threads, $output_dir, @input_files) = @_;

	die if(!defined $bin_fastqc);
	die if(!looks_like_number($threads));
	die if(!defined $output_dir);
	die if(@input_files == 0);

	my %input_files;
	$input_files{$_} = 1 for @input_files;
	die "Error: Ambiguous input files: ".join(", ", @input_files)."\n" if(keys(%input_files) != @input_files);

	my $cmd = "mkdir $output_dir";
	my $cmd2 = "$bin_fastqc -t $threads --outdir $output_dir --extract ".join(" ", @input_files);

	return ($cmd, $cmd2);
}


# -----------------------------------------------
sub _merge_bam_files_cmd {
# -----------------------------------------------
	my ($bin_samtools, $merged_mapping_file, @mapping_files) = @_;

	die "No input files for '_merge_bam_files_cmd' provided\n" if(@mapping_files == 0);

	my $cmd = "$bin_samtools merge $merged_mapping_file ".join(" ", @mapping_files);

	return $cmd;
}


# -----------------------------------------------
sub _picard_metrics_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_picard, $mem, $reference, $input_bam, $output_prefix) = @_;

	die if(!defined $mem);

	my @cmds;
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT R=$reference I=$input_bam O=$output_prefix.alignment_summary.txt &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectBaseDistributionByCycle VALIDATION_STRINGENCY=SILENT I=$input_bam CHART=$output_prefix.base_distr_by_cycle.pdf O=$output_prefix.base_distr_by_cycle_summary.txt &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectGcBiasMetrics VALIDATION_STRINGENCY=SILENT R=$reference I=$input_bam O=$output_prefix.gc_bias_metrics.txt CHART=$output_prefix.gc_bias_metrics.pdf S=$output_prefix.base_distr_by_cycle_summary.txt &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$input_bam O=$output_prefix.insert_size.txt H=$output_prefix.insert_size.pdf M=0.5 &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectQualityYieldMetrics VALIDATION_STRINGENCY=SILENT I=$input_bam O=$output_prefix.quality_yield.txt &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectRawWgsMetrics VALIDATION_STRINGENCY=SILENT R=$reference I=$input_bam O=$output_prefix.raw_wgs.txt INCLUDE_BQ_HISTOGRAM=true &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard CollectWgsMetrics VALIDATION_STRINGENCY=SILENT R=$reference I=$input_bam O=$output_prefix.wgs.txt INCLUDE_BQ_HISTOGRAM=true &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard MeanQualityByCycle VALIDATION_STRINGENCY=SILENT I=$input_bam O=$output_prefix.mean_qual_by_cycle.txt CHART=$output_prefix.mean_qual_by_cycle.pdf &");
	push(@cmds, "$bin_java -Xmx$mem -jar $jar_picard QualityScoreDistribution VALIDATION_STRINGENCY=SILENT I=$input_bam O=$output_prefix.qual_score_dist.txt CHART=$output_prefix.qual_score_dist.pdf &");
	push(@cmds, "wait");

	return @cmds;
}


# -----------------------------------------------
sub _trimmomatic_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_trimmomatic, $adapter_file, $cores, $mem, $fastq1, $fastq2, $output_forward_paired_file, $output_reverse_paired_file, $output_forward_unpaired_file, $output_reverse_unpaired_file, $log_file) = @_;

	die if(!defined $bin_java);
	die if(!defined $jar_trimmomatic || !($jar_trimmomatic =~ /\.jar$/));
	die if(!looks_like_number($cores));
	die if(!defined $mem);
	die if(!defined $output_forward_paired_file);
	die if(!defined $output_reverse_paired_file);
	die if(!defined $output_forward_unpaired_file);
	die if(!defined $output_reverse_unpaired_file);
	die if(!defined $log_file);

	my $gc_cores = ($cores >= 4) ? 4 : $cores;

	# Phred33 or Phred64: https://www.biostars.org/p/63225/
	# Automatic detection in the new version
	# -XX:ParallelGCThreads=$gc_cores
	my $cmd = "$bin_java -XX:ParallelGCThreads=$gc_cores -Xmx$mem -jar $jar_trimmomatic PE -phred33 -threads $cores -trimlog $log_file $fastq1 $fastq2 $output_forward_paired_file $output_forward_unpaired_file $output_reverse_paired_file $output_reverse_unpaired_file LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50";

	if(defined $adapter_file){
		$cmd .= " ILLUMINACLIP:$adapter_file:2:30:10";
	}

	return $cmd;
}


# -----------------------------------------------
sub _verify_bam_id_cmd {
# -----------------------------------------------
	my ($bin_verify_bam_id, $input_vcf, $input_bam, $output_prefix) = @_;

	my $cmd = "$bin_verify_bam_id --vcf $input_vcf --bam $input_bam --out $output_prefix --verbose --best";

	return $cmd;
}


# -----------------------------------------------
sub _mapping_cmd {
# -----------------------------------------------
	my ($read_group, $bwa_bin, $bin_samtools, $bwa_index, $fastq_r1, $fastq_r2, $output_file, $n_cores, $mem_samtools) = @_;

	die if(!defined $read_group || $read_group eq "");
	die if(!defined $bwa_bin);
	die if(!defined $bwa_index);
	die if(!defined $fastq_r1);
	die if(!defined $fastq_r2);
	die if (!($output_file =~ /bam$/));
	die if (!looks_like_number($n_cores));

	# Don't use -@ or -mem when piping out of samtools view

	my $cmd = "$bwa_bin mem -R '".join("\\t", ('@RG', $read_group))."' -M -t $n_cores $bwa_index $fastq_r1 $fastq_r2 | $bin_samtools sort -@ $n_cores -m $mem_samtools -o $output_file -";

	return $cmd;
}


# -----------------------------------------------
sub _mark_duplicates {
# -----------------------------------------------
	my ($bin_java, $mem, $cores, $jar_picard, $input_file, $output_file, $duplication_metrics_file, $tmp_dir) = @_;

	my $cmd = "$bin_java -XX:ParallelGCThreads=$cores -Xmx$mem -jar $jar_picard MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=$input_file OUTPUT=$output_file M=$duplication_metrics_file TMP_DIR=$tmp_dir";

	return $cmd;
}


# -----------------------------------------------
sub _base_recal_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $cores, $mem, $reference, $known_sites, $input_bam, $recal_report_file) = @_;

	die if(!looks_like_number($cores));
	die if(!defined $mem);
	die if(keys(%{$known_sites}) == 0);

	#https://software.broadinstitute.org/gatk/documentation/article?id=2801
	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk -T BaseRecalibrator -nct $cores -R $reference -I $input_bam -ics 8 -mcs 4 -o $recal_report_file -knownSites ".join(" -knownSites ", keys(%{$known_sites}));

	return ($cmd);
}


# -----------------------------------------------
sub _apply_base_recal_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $cores, $mem, $reference, $input_bam, $recal_report_file, $output_bam) = @_;

	die if(!looks_like_number($cores));
	die if(!defined $mem);

	#https://software.broadinstitute.org/gatk/documentation/article?id=2801
	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk -T PrintReads -nct $cores -R $reference --filter_mismatching_base_and_quals -I $input_bam -BQSR $recal_report_file -o $output_bam";

	return ($cmd);
}


# -----------------------------------------------
sub _base_recal_stats_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $cores, $mem, $reference, $known_sites, $input_bam, $recal_report_file, $post_recal_report_file, $recal_plots) = @_;

	die if(!looks_like_number($cores));
	die if(!defined $mem);
	die if(keys(%{$known_sites}) == 0);

	#https://software.broadinstitute.org/gatk/documentation/article?id=2801
	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk -T BaseRecalibrator -nct $cores -R $reference -I $input_bam -ics 8 -mcs 4 -BQSR $recal_report_file -o $post_recal_report_file -knownSites ".join(" -knownSites ", keys(%{$known_sites}));
	my $cmd2 = "$bin_java -Xmx$mem -jar $jar_gatk -T AnalyzeCovariates -R $reference -before $recal_report_file -after $post_recal_report_file -plots $recal_plots";

	return ($cmd, $cmd2);
}


# -----------------------------------------------
sub _bam_index_cmd {
# -----------------------------------------------
	my ($bin_samtools, $file) = @_;

	return "$bin_samtools index $file";
}


# -----------------------------------------------
sub _vcf_index_cmd {
	my ($bin_tabix, $file) = @_;

	return "$bin_tabix -f -p vcf $file.g.vcf.gz";
}


# -----------------------------------------------
sub _hc_per_sample_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $reference, $interval_list, $cores, $mem, $input_file, $output_file) = @_;

	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/) && !($reference =~ /fa.gz$/) && !($reference =~ /fasta.gz$/));
	die if(!looks_like_number($cores));
	die if (!($jar_gatk =~ /jar$/));
	die $input_file if (!($input_file =~ /bam$/));
	die if (!($output_file =~ /g\.vcf\.gz$/));


	my $cmd = "$bin_java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 ";
	$cmd .= (defined $mem && $mem ne "") ? "-Xmx$mem " : "";
	$cmd .= "-jar $jar_gatk \\
	-R $reference \\
	-T HaplotypeCaller -nct $cores --emitRefConfidence BP_RESOLUTION \\
	--variant_index_type LINEAR --variant_index_parameter 128000 \\
	-I $input_file \\
	-o $output_file";


	if(defined $interval_list && $interval_list ne ""){
		$cmd .= " -L $interval_list";
	}


	return $cmd;
}


## -----------------------------------------------
#sub _merge_gvcfs_cmd {
## -----------------------------------------------
#	my ($bin_java, $jar_gatk, $bin_bgzip, $reference, $cores, $mem, $input_file, $output_file) = @_;
#
#	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
#	die if (!($jar_gatk =~ /jar$/));
#	die if (!-e $input_file);
#	die if (!($output_file =~ /g\.vcf\.gz$/));
#	die if(!looks_like_number($cores));
#	die if(!defined $mem);
#
#	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk \\
#	-R $reference \\
#	-T CombineGVCFs \\
#	--variant $input_file \\
#	-o /dev/stdout | $bin_bgzip -c >$output_file";
#
#	return $cmd;
#}


# -----------------------------------------------
sub _parallel_merge_gvcfs_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $parallelize_cmds, $reference, $chromosomes, $cores, $mem, $input_files, $output_file_pref) = @_;


	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
	die if (!($jar_gatk =~ /jar$/));
	die if(!defined $mem);

	my @all_cmds;

	my @select_cmds;
	my @combine_cmds;
	my @tmp_files2;
	foreach my $chr (@{$chromosomes}){

		my @tmp_files;

		for(my $i=0; $i<@{$input_files}; $i++){

			my $input_file = ${$input_files}[$i];
			my $tmp_file = "$output_file_pref.$chr.$i.g.vcf.gz";
			push(@tmp_files, $tmp_file);

			push(@select_cmds, "$bin_java -Xmx$mem -jar $jar_gatk -T SelectVariants -nt 1 -R $reference -V $input_file -o $tmp_file -L $chr")

		}

		my $tmp_file = ${output_file_pref}.".$chr.g.vcf.gz";
		push(@tmp_files2, $tmp_file);

		push(@combine_cmds, "$bin_java -Xmx$mem -jar $jar_gatk -R $reference -T CombineGVCFs -V ".join(" -V ", @tmp_files)." -o $tmp_file");

	}


	push(@all_cmds, "$parallelize_cmds $cores '".join("' \\\n'", @select_cmds)."'\n");
	push(@all_cmds, "$parallelize_cmds $cores '".join("' \\\n'", @combine_cmds)."'\n");
	push(@all_cmds, "$bin_java -Xmx$mem -cp $jar_gatk org.broadinstitute.gatk.tools.CatVariants -R $reference --assumeSorted -V ".join(" -V ", @tmp_files2)." -out ${output_file_pref}.g.vcf.gz\n");


	return @all_cmds;
}


# -----------------------------------------------
sub _cov_depth_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $cores, $mem, $reference, $gene_list, $input_bam_list, $output_suffix) = @_;

	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
	die if (!($jar_gatk =~ /jar$/));
	die if(!looks_like_number($cores));
	die if(!defined $mem);

	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk -T DepthOfCoverage \\
				-R $reference -I $input_bam_list -o $output_suffix \\
				-ct 4 -ct 6 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -ct 35 -ct 40 \\
				-pt sample -pt readgroup -pt library \\
				-geneList $gene_list";

	return ($cmd);
}


# -----------------------------------------------
sub _multi_sample_hc_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $reference, $dbsnp_file, $cores, $mem, $input_files, $output_file) = @_;

	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
	die if (!($jar_gatk =~ /jar$/));
	die if (!($output_file =~ /vcf\.gz$/));
	die if(!looks_like_number($cores));
	die if (!($dbsnp_file =~ /vcf\.gz$/));
	die if(!defined $mem);

	my $gc_cores = ($cores >= 4) ? 4 : $cores;


	#TODO: -nt produces error
	my $cmd = "$bin_java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=$gc_cores ";
	$cmd .= (defined $mem && $mem ne "") ? "-Xmx$mem " : "";
	$cmd .= "-jar $jar_gatk \\
	-nt 1 \\
	--dbsnp $dbsnp_file \\
	-R $reference \\
	-T GenotypeGVCFs \\
	-o $output_file";


	$cmd = $cmd." --variant ".join(" --variant ", @{$input_files});

	return $cmd;
}


# -----------------------------------------------
sub _recal_snv_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $cores, $mem, $reference, $resources, $input_file, $outfile_recal, $outfile_tranches, $rscript_file, $output_file) = @_;

	die if (!($jar_gatk =~ /jar$/));
	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
	die if (!($input_file =~ /vcf\.gz$/ || $input_file =~ /vcf$/));
	die if (!($output_file =~ /vcf\.gz$/ || $output_file =~ /vcf$/));
	die if(!looks_like_number($cores));
	die if(!defined $mem);

	my $gc_cores = ($cores >= 4) ? 4 : $cores;

	#TODO: -nt $cores gives error
	my $cmd = "$bin_java -XX:ParallelGCThreads=$gc_cores -Xmx$mem -jar $jar_gatk \\
				-R $reference \\
				-T VariantRecalibrator \\
				-input $input_file \\
				-recalFile $outfile_recal -tranchesFile $outfile_tranches -rscriptFile $rscript_file\\
				-nt 1 \\
				-mode SNP \\
				-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP \\
				--maxGaussians 6";


				foreach my $res (keys(%{$resources})){

					my $known = ${$resources}{$res}{'known'}; die if (!defined $known || $known eq "");
					my $training = ${$resources}{$res}{'training'}; die if (!defined $training || $training eq "");
					my $truth = ${$resources}{$res}{'truth'}; die if (!defined $truth || $truth eq "");
					my $prior = ${$resources}{$res}{'prior'}; die if (!defined $prior || $prior eq "");
					die "Prior '$prior' not a number" if(!looks_like_number($prior));
					my $file = ${$resources}{$res}{'file'};

					$cmd .= " -resource:$res,known=$known,training=$training,truth=$truth,prior=$prior $file";
				}


	my $cmd2 = "$bin_java -XX:ParallelGCThreads=$gc_cores -Xmx$mem -jar $jar_gatk \\
				-R $reference \\
				-T ApplyRecalibration \\
				-input $input_file -recalFile $outfile_recal -tranchesFile $outfile_tranches \\
				-o $output_file \\
				-mode SNP --ts_filter_level 99.0";


	return ($cmd, $cmd2);
}


# -----------------------------------------------
sub _recal_indel_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $cores, $mem, $reference, $resources, $input_file, $outfile_recal, $outfile_tranches, $rscript_file, $output_file) = @_;

	die if (!($jar_gatk =~ /jar$/));
	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
	die if (!($input_file =~ /vcf\.gz$/ || $input_file =~ /vcf$/));
	die if (!($output_file =~ /vcf\.gz$/ || $output_file =~ /vcf$/));
	die if(!looks_like_number($cores));
	die if(!defined $mem);

	my $gc_cores = ($cores >= 4) ? 4 : $cores;


	#TODO: -nt $cores gives error
	my $cmd = "$bin_java -XX:ParallelGCThreads=$gc_cores -Xmx$mem -jar $jar_gatk \\
				-R $reference \\
				-T VariantRecalibrator \\
				-input $input_file \\
				-recalFile $outfile_recal -tranchesFile $outfile_tranches -rscriptFile $rscript_file\\
				-mode INDEL \\
				-nt 1 \\
				-an FS -an ReadPosRankSum -an MQRankSum -an DP \\
				--maxGaussians 4 --minNumBadVariants 5000";

				foreach my $res (keys(%{$resources})){

					my $known = ${$resources}{$res}{'known'}; die if (!defined $known || $known eq "");
					my $training = ${$resources}{$res}{'training'}; die if (!defined $training || $training eq "");
					my $truth = ${$resources}{$res}{'truth'}; die if (!defined $truth || $truth eq "");
					my $prior = ${$resources}{$res}{'prior'}; die if (!defined $prior || $prior eq "");
					die "Prior '$prior' not a number" if(!looks_like_number($prior));
					my $file = ${$resources}{$res}{'file'};

					$cmd .= " -resource:$res,known=$known,training=$training,truth=$truth,prior=$prior $file";
				}

	my $cmd2 = "$bin_java -XX:ParallelGCThreads=$gc_cores -Xmx$mem -jar $jar_gatk \\
				-R $reference \\
				-T ApplyRecalibration \\
				-input $input_file -recalFile $outfile_recal -tranchesFile $outfile_tranches \\
				-o $output_file \\
				-mode INDEL --ts_filter_level 99.0";

	return ($cmd, $cmd2);
}


# -----------------------------------------------
sub _calc_gen_post_cmd {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $mem, $reference, $gold_genotypes, $input_file, $output_file_cgp, $output_file) = @_;

	die if (!($jar_gatk =~ /jar$/));
	die if (!($reference =~ /fa$/) && !($reference =~ /fasta$/));
	die if (!($input_file =~ /vcf\.gz$/ || $input_file =~ /vcf$/));
	die if (!($output_file =~ /vcf\.gz$/ || $output_file =~ /vcf$/));
	die if(!defined $mem);

	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk -T CalculateGenotypePosteriors \\
				-R $reference  --supporting $gold_genotypes -V $input_file -o $output_file_cgp";


	my $cmd2 = "$bin_java -Xmx$mem -jar $jar_gatk -T VariantFiltration \\
				-R $reference -V $output_file_cgp -G_filter \"GQ < 20.0\" -G_filterName lowGQ -o $output_file";

	return ($cmd, $cmd2);
}


# -----------------------------------------------
sub _create_clean_variant_sets {
# -----------------------------------------------
	my ($bin_java, $jar_gatk, $mem, $reference, $input_file, $output_vars, $output_snv, $output_indel) = @_;


	#$pl_filterVCFbyAD -g 30 -d 10 -o $output_file -rn
	my $cmd = "$bin_java -Xmx$mem -jar $jar_gatk -T SelectVariants -R $reference -V $input_file --excludeFiltered -o $output_vars";
	my $cmd2 = "$bin_java -Xmx$mem -jar $jar_gatk -T SelectVariants -R $reference -V $input_file --selectTypeToExclude INDEL --excludeFiltered -o $output_snv";
	my $cmd3 = "$bin_java -Xmx$mem -jar $jar_gatk -T SelectVariants -R $reference -V $input_file -selectType INDEL --excludeFiltered -o $output_indel";

	return ($cmd, $cmd2, $cmd3);
}


# -----------------------------------------------
sub _file_exists {
# -----------------------------------------------
	my ($size, @output_files) = @_;

	foreach my $output_file (@output_files){
		if (-e $output_file && -d $output_file){

			my $dir = dir($output_file);
			my $is_empty = $dir -> stat && !$dir -> children;

			if($is_empty){
				return 0;
			}
		}
		elsif (!-e $output_file || -s $output_file <= $size){
			return 0;
		}
	}

	return 1;
}


# -----------------------------------------------
sub _is_successful {
# -----------------------------------------------
	my ($log) = @_;


	if (-e $log){

		local $/ = undef;
		open(IN, "<".$log) or die "Cannot open file '$log': $!\n";
		my $log_string = <IN>;
		close IN;


		if ($log_string =~ /ERROR/i){
			return (0, "ERROR");
		}
		elsif ($log_string =~ /EXCEPTION/i){
			return (0, "EXCEPTION");
		}
		elsif ($log_string =~ /CANCELLED/i){
			return (0, "CANCELLED");
		}
		elsif($log_string =~ /No space left on device/i){
			return (0, "No space left on device");
		}
		elsif($log_string =~ /Disk quota exceeded/i){
			return (0, "Disk quota exceeded");
		}
		elsif($log_string =~ /slurmstepd: Exceeded job memory limit/i){
			return (0, "slurmstepd: Exceeded job memory limit");
		}
		elsif($log_string =~ /exited with non-zero status 1/i){
			return (0, "exited with non-zero status 1");
		}
		elsif($log_string =~ /\@COMPLETE$/i) {
			return (1, "");
		}
		else {
			return (0, "unknown");
		}
	}
	else {
		return (0, "log file doesn't exist");
	}
}


# -----------------------------------------------
sub _write_to_file {
# -----------------------------------------------
	my ($string, $output_file) = @_;

	die if (!defined $string);
	die if (!-e dirname($output_file));

	open(OUT, ">".$output_file) or die "Cannot open file '$output_file': $!\n";
	print OUT $string;
	close OUT;

	return 1;
}


# -----------------------------------------------
sub _time {
# -----------------------------------------------
	my (@cmds) = @_;

	my @time_cmds;
	foreach my $c (@cmds){

		if(defined $c && $c eq "wait"){
			push(@time_cmds, $c);
		}
		elsif(defined $c && $c =~ /^cd /){
			push(@time_cmds, $c);
		}
		elsif(defined $c && $c =~ /^echo /){
			push(@time_cmds, $c);
		}
		elsif(defined $c && $c ne ""){
			push(@time_cmds, "/usr/bin/time -v $c");
		}
		else {
			die;
		}
	}

	return @time_cmds;
}


# -----------------------------------------------
sub _path {
# -----------------------------------------------
	my ($path, @files) = @_;

	my @results;
	foreach (@files){
		push(@results, $path."/$_") if (defined $_);
	}

	return @results;
}


# -----------------------------------------------
sub _submit_slurm_job {
# -----------------------------------------------
	my ($cmds, $slurm_file, $slurm_pars) =  @_;

	my $slurm_script = _create_job_script($cmds, $slurm_pars);
	_write_to_file($slurm_script, $slurm_file) or die;

	my $output = `sbatch $slurm_file`;

	if ($output =~ /Submitted batch job [0-9]+ /){
		$output =~ /Submitted batch job ([0-9]+) /;
		my $job_id = $1;
		die if (!looks_like_number($job_id));
		return $job_id;
	}
	else {
		die "Error occured when submitting to slurm: ".$output."\n";
	}
}


# -----------------------------------------------
sub _cp_cmd {
# -----------------------------------------------
	my $target = pop @_;
	my @source = @_;

	my @cmds;
	foreach my $s (@source){
		push(@cmds, "cp -r $s $target");
	}

	return @cmds;
}


# -----------------------------------------------
sub _create_output_dir {
# -----------------------------------------------
	my ($ids, $output_root_dir) = @_;

	my $output_dir = $output_root_dir;
	$output_dir .= "/".join("/", @{$ids}) if(@{$ids} > 0); system("mkdir -p $output_dir") if (!-e $output_dir);

	return ($output_dir);
}


# -----------------------------------------------
sub _create_job_script {
# -----------------------------------------------
	my ($cmds, $slurm_pars) = @_;

	my @slurm_pars;
	foreach my $par (sort keys(%{$slurm_pars})){

		if (defined ${$slurm_pars}{$par}){
			push(@slurm_pars, "#SBATCH $par=".${$slurm_pars}{$par});
		}
		else {
			push(@slurm_pars, "#SBATCH $par");
		}
	}
	my $slurm_pars_string = join("\n", @slurm_pars);

	my $cmds_string = join("\n", @{$cmds});

	my $script = <<"END_SCRIPT";
#! /bin/bash
$slurm_pars_string

$cmds_string

END_SCRIPT

	return $script;
}


# -----------------------------------------------
sub _run {
# -----------------------------------------------
	my ($cmds, $slurm_file, $log_file, $output_files, $continue, $cpus, $mem, $time, $dependencies, $mode) = @_;


	die "CPUs: $cpus\n" if(!looks_like_number($cpus));
	die "Mem: $mem\n" if(!($mem =~ /GB$/));
	die "Time: $time\n" if(!looks_like_number($time));
	die if(ref($dependencies) ne 'ARRAY');
	die if(!looks_like_number($continue));
	die if(!-e dirname($slurm_file));


	my ($is_succ, $message) = _is_successful($log_file);


	if($mode eq "error"){

		if (!$is_succ){
			print STDERR "Error(s) in '$log_file'\n";
			print STDERR "Reason: $message\n";
		}
		elsif(!_file_exists(1, @{$output_files})){
			print STDERR "One or more files do not exist: ".join(", ", @{$output_files})."\n";
		}

		return (undef, $continue);
	}


	if ($continue == 1 && !$is_succ){
		print STDERR "Error(s) in '$log_file'. Re-run job\n" if(-e $log_file);
	}
	elsif($continue == 1 && !_file_exists(1, @{$output_files})){
		print STDERR "File(s) ".join(",", @{$output_files})." do not exist or are empty. Re-run job\n";
	}
	elsif($continue == 1) {
		return (undef, $continue);
	}


	my @nodes_excl;
	my @partitions;
	if($time <= 14*24*60){
		push(@partitions, "longterm");
		#push(@partitions, "lied");
	}
	if($time <= 3*24*60){
		push(@partitions, "shortterm");
	}
	if($time <= 12*60){
		push(@partitions, "debug");
	}
	if($time > 14*24*60){
		die;
	}


	$time = _convert_minutes_to_slurm_time($time);


	my %slurm_pars = ("--cluster" => "omics",
						"--partition" => join(",", reverse @partitions),
						"--exclude" => join(",", @nodes_excl),
						"--requeue" => undef,
						"--cpus-per-task" => $cpus,
						"--mem" => $mem,
						"--time" => $time,
						"--output" => $log_file);
	if (@{$dependencies} > 0){

		my @d;
		foreach my $d (@{$dependencies}){
			if(looks_like_number($d)){
				push(@d, $d);
			}
		}
		$slurm_pars{"--dependency"} = "afterok:".join(":", @d) if(@d > 0);
		$slurm_pars{"--kill-on-invalid-dep"} = "yes";
	}


	if(-e $log_file){
		unlink($log_file) or die "An error occurred when deleting '$log_file'\n";
	}


	my $job_id = _submit_slurm_job([_time(@{$cmds}), 'echo @COMPLETE'], $slurm_file, \%slurm_pars);
	print $job_id.": ".$slurm_file."\n";


	return ($job_id, 0);
}


# -----------------------------------------------
sub _convert_minutes_to_slurm_time {
# -----------------------------------------------
	my ($total_min) = @_;

	my $days = int($total_min/(24*60));
	my $hours = ($total_min/60)%24;
	my $min = $total_min%60;

	return "$days-$hours:$min:0";
}


# -----------------------------------------------
sub _basename_list {
# -----------------------------------------------
	my @l = @_;

	my @results;
	push(@results, basename($_)) for @l;

	return @results;
}


# -----------------------------------------------
sub _get_best_split {
# -----------------------------------------------
	my ($arr, $min, $max) = @_;

	die "Error: Wrong parameters set (min=$min and max=$max)\n" if($min > $max);

	my $l = @{$arr};
	if($l < $min){
		print STDERR "Warning: Length of array ($l) less than minimum allowed split size (min=$min). Skip splitting\n";
		return ($arr);
	}
	elsif($max > $l) {
		$max = $l;
	}

	my $n_splits;
	my $orig_max = $max;
	while(){

		if($max < $min) {
			die "Error: can't split with parameters set (min=$min and max=$orig_max)\n";
		}
		elsif($l % $max >= $min){

			$n_splits = floor($l / $max) + 1;
			last;
		}
		elsif($l % $max == 0){

			$n_splits = floor($l / $max);
			last;
		}
		else {
			$max = $max-1;
		}
	}


	my @results;
	my $i = 0;
	foreach my $e (@{$arr}){

		if ($i < $n_splits){
			push(@{$results[$i]}, $e);
		}
		else {
			$i = 0;
			push(@{$results[$i]}, $e);
		}

		$i++;
	}

	return @results;
}
