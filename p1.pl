use strict;
use warnings;
use POSIX;

########################This is the script for discovering variants based on RNA-seq data########################
############The input file should be an RNA-seq bam file generated using the STAR pipeline (/STAR.sh)############
############The input file should be an RNA-seq bam file generated using the STAR pipeline (/STAR.sh)############

my $start_time_epoch = time; #epoch time 
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $start_time_str = sprintf("%02d-%02d-%02d_%02d-%02d-%02d", $year-100,$mon+1, $mday, $hour, $min, $sec);

my $input = join "\t", @ARGV;
$input =~ s/\s+$//;
unless ((
		#input_files
		$input =~ /-Tumor\t/
		&& $input =~ /-Normal\t/
		&& $input =~ /-R\t/
		&& $input =~ /-outPrefix\t/
		
		#tools
		&& $input =~ /-gatk\t/
		&& $input =~ /-picard\t/
		&& $input =~ /-hisat2\t/
		
		#reference
		&& $input =~ /-gtf\t/
		&& $input =~ /-dbsnp\t/
		&& $input =~ /-hisat2_reference\t/
		)		
		|| $input eq '-h'){
	print "Found missing paramters, Please use -h for help information\n";
	exit;	
}

if ($input eq '-h'){
	print "This is the script for discovering variants based on RNA-seq data. The input file should be an RNA-seq bam file generated using the STAR pipeline.\n";
	print "Usage: perl [options]... -mode single -d ./data_files -o ./out_files -g ./GRCH38 -f pair-end\n";
	print "The result will be stored in output directory, name with Final_*";
	print "###version: 1.0.0\n";
	print "#########################################################################################################################################\n";
	print "Required:\n";	
	print "-I\n\t./data_files: The directory for input file\n\n";
	print "-O\n\t./out_files: The directory for output file\n\n";
	print "-R\n\t./genome.fasta: Reference sequence file\n\n";
	
	exit;
}

###############################################################################
#prepare the input data from arugument input
my $tumorInput = '.';
if ($input =~ /-Tumor\t(\S+)/){
	$tumorInput = $1;
}

my $normalInput = '.';
if ($input =~ /-Normal\t(\S+)/){
	$normalInput = $1;
}

my $fastaReference = '.';
if ($input =~ /-R\t(\S+)/){
	$fastaReference = $1;
}

my $outDir = '.';
if ($input =~ /-outPrefix\t(\S+)/){
	$outDir = $1;
}

###############################################################################
my $gtfFile = '.';
if ($input =~ /-gtf\t(\S+)/){
	$gtfFile = $1;
}

my $dbsnpFile = '.';
if ($input =~ /-dbsnp\t(\S+)/){
	$dbsnpFile = $1;
}

my $hisat2Reference = '.';
if ($input =~ /-hisat2_reference\t(\S+)/){
	$hisat2Reference = $1;
}

###############################################################################
my $gatk = '.';
if ($input =~ /-gatk\t(\S+)/){
	$gatk = $1;
}

my $picard = '.';
if ($input =~ /-picard\t(\S+)/){
	$picard = $1;
}

my $hisat2 = '.';
if ($input =~ /-hisat2\t(\S+)/){
	$hisat2 = $1;
}

###########################################Generate output directory#############################################
unless(-d $outDir){
	system("mkdir $outDir");
}

##############################################Preprocess bam file################################################
my $ram = 50;
RNAPreProcessParallel($outDir, $tumorInput);


my %gtf = readGTF($gtfFile);



###################################################Subroutines###################################################
sub readGTF{
	my ($filename) = @_;
	my %gtf;	
	if(-e $filename){
		open(GTF, "$filename") or die "Cannot open $filename for reading: $!\n";
		while (<GTF>) {
			next if $_ =~ /^\#/;
			$_ =~ s/\s+$//;
			my @line = split /\t/, $_;
			next unless $line[2] eq 'gene';
			$gtf{$line[0]}{$line[3]} = 1;
			$gtf{$line[0]}{$line[4]} = 1;
		}
		close(GTF);
		unless (keys %gtf) {
			print "Can't find the gene information in $filename, please check your gtf annotation file\n";
			exit;
		} 
	}
	else{
		print "Can't find the GTF annotation file in $filename\n";
		exit;
	}
	return  %gtf;
}

sub RNAPreProcessParallel{
	my ($outDir,$bamFile) = @_;
	
	print "##################################################################################################################################################\n";	
	system("rm $outDir/GATK.log.out") if -e "$outDir/GATK.log.out";
	system("rm -r $outDir/tmp") if -d "$outDir/tmp";
	my $start_time_epoch = time; #epoch time 	
	my $start_time_str = getTime();
	
	print "$start_time_str START Pre-process for $bamFile\n";
	system("mkdir $outDir/tmp") unless -d "$outDir/tmp";
	system("java -jar $picard SplitSamByNumberOfReads -I $bamFile -O $outDir/tmp -N_READS 10000000 -VALIDATION_STRINGENCY LENIENT 2>> $outDir/GATK.log.out");
	my $upLimit = floor($ram/20);	
	my $fileCount = `ls $outDir/tmp/shard_*.bam | wc -l`;
	$fileCount =~ s/\s+$//;	
	my $iteration = ceil($fileCount/$upLimit);
	foreach my $i(1..$iteration){
		my $markDupCommand = "";
		foreach my $j(1..$upLimit){
			my $temp = ($i-1)*$upLimit+$j;
			last if $temp > $fileCount;
			$temp = sprintf("%04d", $temp);
			$markDupCommand = $markDupCommand . "java -jar $picard MarkDuplicates -I $outDir/tmp/shard_$temp.bam -O $outDir/tmp/md_shard_$temp.bam -M $outDir/tmp/md_shard_$temp.bam_md_metrics.txt --REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT 2>> $outDir/GATK.log.out &\n";
		}
		$markDupCommand = $markDupCommand . "wait\n";		
		system($markDupCommand);
		
	}
	
	$upLimit = floor($ram/10);	
	$iteration = ceil($fileCount/$upLimit);
	foreach my $i(1..$iteration){
		my $splitCommand = "";
		foreach my $j(1..$upLimit){
			my $temp = ($i-1)*$upLimit+$j;
			last if $temp > $fileCount;
			$temp = sprintf("%04d", $temp);
			$splitCommand = $splitCommand . "$gatk SplitNCigarReads -R $fastaReference -I $outDir/tmp/md_shard_$temp.bam -O $outDir/tmp/split_md_shard_$temp.bam 2>> $outDir/GATK.log.out &\n";
		}
		$splitCommand = $splitCommand . "wait\n";		
		system($splitCommand);
		
	}
	my $BRCommand = "$gatk BaseRecalibrator -R $fastaReference ";
	foreach my $i(1..$fileCount){
		$i = sprintf("%04d", $i);
		$BRCommand = $BRCommand . "-I $outDir/tmp/split_md_shard_$i.bam ";		
	}
	$BRCommand = $BRCommand . "-known-sites $dbsnpFile -O $outDir/tmp/Aligned.out.table 2>> $outDir/GATK.log.out";	
	system($BRCommand);	
	
	foreach my $i(1..$iteration){
		my $BQSRCommand = "";
		foreach my $j(1..$upLimit){
			my $temp = ($i-1)*$upLimit+$j;
			last if $temp > $fileCount;
			$temp = sprintf("%04d", $temp);
			$BQSRCommand = $BQSRCommand . "$gatk ApplyBQSR -R $fastaReference -I $outDir/tmp/split_md_shard_$temp.bam -bqsr-recal-file $outDir/tmp/Aligned.out.table -O $outDir/tmp/re_split_md_shard_$temp.bam 2>> $outDir/GATK.log.out &\n";
		}
		$BQSRCommand = $BQSRCommand . "wait\n";		
		system($BQSRCommand);
		
	}
	
	my $mergeCommand = "java -jar $picard MergeSamFiles ";
	foreach my $i(1..$fileCount){
		$i = sprintf("%04d", $i);
		$mergeCommand = $mergeCommand . "-I $outDir/tmp/re_split_md_shard_$i.bam ";		
	}
	$mergeCommand = $mergeCommand . "-O $outDir/reAligned.bam 2>> $outDir/GATK.log.out";	
	system($mergeCommand);
	system("rm -r $outDir/tmp");
	system("samtools index $outDir/reAligned.bam");
	
	my $end_time_str = getTime();
	my $end_time_epoch = time;
	my $pre_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str END Pre-process for $bamFile\n";
	print "$end_time_str Time Period: ", $pre_epoch, "\n";
	print "##################################################################################################################################################\n";
}

sub variantsCalling{
	my ($tumorInput,$normalInput,$ID,$outDir,$genomeDir) = @_;
	my $start_time_epoch = time; #epoch time
	my $start_time_str = getTime();
	print "##################################################################################################################################################\n";
	print "$start_time_str START $ID Variants calling\n";
	
	#test if bed file exist
	unless (-e $gtfBedDir){
		system("awk '{print \$1:\$4-\$5}' $gtfBedDir | grep -v '#' > $genomeDir/GRCH38_bed.list");
		$gtfBedDir = "$genomeDir/GRCH38_bed.list";
	}
	my $lineCount = `wc -l $gtfBedDir`;
	$lineCount = $1 if $lineCount =~ /(\d+)\s/;
	my $lineBreaker = ceil($lineCount/$thread);
	
	system('rm -r $genomeDir/bedSplit*.list');
	system("split -l $lineBreaker $gtfBedDir $genomeDir/bedSplit");
	my @bedList = `ls $genomeDir/bedSplit*`;
	my $command = '';	
	
	my $tumor_ID_string = `samtools view -H $tumorInput | grep SM: | head -1`;
	my $tumor_ID = $1 if $tumor_ID_string =~ /SM\:(\S+)\s/;
	print "$tumor_ID\n";
	
	my $normal_ID_string = `samtools view -H $normalInput | grep SM: | head -1`;
	my $normal_ID = $1 if $normal_ID_string =~ /SM\:(\S+)\s/;
	print "$normal_ID\n";
	
	system("rm $outDir/$ID\_*.log.out");
	foreach my $bedList (@bedList){
		$bedList =~ s/\s+$//;
		next if $bedList =~ /vcf/;
		next if $bedList =~ /list/;
		system("mv $bedList $bedList.list");
		my $bed = $1 if $bedList =~ /$genomeDir\/(\S+)/;
		print "$bedList\n";
		print "$bed\n";			
		$command = $command . "$gatk Mutect2 -R $fastaDir -I $tumorInput -I $normalInput -tumor $tumor_ID -normal $normal_ID -O $outDir/$ID\_$bed.vcf -L $bedList.list --panel-of-normals /projects/pharma_xwang317/mcho43/genome/somatic-hg38_1000g_pon.hg38.vcf.gz --germline-resource /projects/pharma_xwang317/mcho43/genome/somatic-hg38_af-only-gnomad.hg38.vcf.gz --read-validation-stringency LENIENT --disable-read-filter AllowAllReadsReadFilter 2>> $outDir/$ID\_$bed.log.out &\n";
	}
	$command = $command . "wait\n";
	print "$command\n";	
	system($command);
	
	my %vcf;
	my @vcf = `ls $outDir/$ID*.vcf`;
	foreach my $vcf (@vcf){
		$vcf =~ s/\s+$//;
		print "$vcf\n";
		open(VCF, "$vcf") or die "Cannot open file for reading: $!\n";
		while (<VCF>) {
			next if $_ =~ /^\#/;	
			$_ =~ s/\s+$//;
			my @line = split /\t/, $_;
			if ($gtfFlag==1){
				if(length $line[3] > 5 or length $line[4] > 5){
					if(exists $gtf{$line[0]}{$line[1]} 
					or exists $gtf{$line[0]}{$line[1]+1} 
					or exists $gtf{$line[0]}{$line[1]+2} 
					or exists $gtf{$line[0]}{$line[1]+3} 
					or exists $gtf{$line[0]}{$line[1]+4} 
					or exists $gtf{$line[0]}{$line[1]+5} 
					or exists $gtf{$line[0]}{$line[1]-1} 
					or exists $gtf{$line[0]}{$line[1]-2} 
					or exists $gtf{$line[0]}{$line[1]-3} 
					or exists $gtf{$line[0]}{$line[1]-4} 
					or exists $gtf{$line[0]}{$line[1]-5} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])+1} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])+2} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])+3} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])+4} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])+5} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])-1} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])-2} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])-3} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])-4} 
					or exists $gtf{$line[0]}{$line[1]+length($line[3])-5} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])+1} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])+2} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])+3} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])+4} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])+5} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])-1} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])-2} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])-3} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])-4} 
					or exists $gtf{$line[0]}{$line[1]+length($line[4])-5} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])+1} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])+2} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])+3} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])+4} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])+5} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])-1} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])-2} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])-3} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])-4} 
					or exists $gtf{$line[0]}{$line[1]-length($line[3])-5} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])+1} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])+2} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])+3} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])+4} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])+5} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])-1} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])-2} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])-3}
					or exists $gtf{$line[0]}{$line[1]-length($line[4])-4} 
					or exists $gtf{$line[0]}{$line[1]-length($line[4])-5})
					{}
					else{
						$vcf{$line[0]}{$line[1]} = $_;
					}
				}
				else{
					$vcf{$line[0]}{$line[1]} = $_;
				}
			}
			else{
				$vcf{$line[0]}{$line[1]} = $_;
			}
		}
		close(VCF);
	}

	my $variantFile = "$outDir/$ID\_Variants.vcf";
	my $variantBed = "$outDir/$ID\_Variants.bed";
	open(OUTFILE, ">$variantFile") or die "Cannot open file for reading: $!\n";
	open(OUTBED, ">$variantBed") or die "Cannot open file for reading: $!\n";
	foreach my $chr(sort {lc $a cmp lc $b} keys %vcf){
		foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
			print OUTFILE "$vcf{$chr}{$pos}\n";
			print OUTBED "$chr\t$pos\t$pos\n";
		}
	}
	close(OUTFILE);
	close(OUTBED);
	
	my $end_time_epoch = time;
	my $end_time_str = getTime();
	my $first_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str END $ID Variants calling\n";
	print "$end_time_str Time Period: ", $first_epoch, "\n";
	print "##################################################################################################################################################\n";
}

sub extractReadFromBam{
	my ($dir, $variantsFile, $bamInput, $tumorSM, $tumorID, $tumorPL, $tumorLB, $tumorPU) = @_;
	my $start_time_epoch = time; #epoch time
	my $start_time_str = getTime();
	print "##################################################################################################################################################\n";
	print "$start_time_str START Extract Reads from $dir Variants Calling\n";

	system("samtools view -L $variantsFile -@ $thread $bamInput | cut -f1 > $dir/IDs_all.txt");
	system("java -Xmx7g -jar $picard FilterSamReads I=$bamInput O=$dir/tmp_bam.bam READ_LIST_FILE=$dir/IDs_all.txt FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=LENIENT 2>> $dir/HISAT2.log.out");
	system("rm $dir/IDs_all.txt");
	system("rm $bamInput");
	system("samtools view -H $dir/tmp_bam.bam | sed '\$d' - > $dir/tmp_header_T.sam");
	system("samtools view $dir/tmp_bam.bam | awk '\$2 < 2040 { print }' > $dir/tmp0_T.sam");
	system("rm $dir/tmp_bam.bam");
	system("cat $dir/tmp_header_T.sam $dir/tmp0_T.sam > $dir/tmp_filteredbamT.sam");
	system("rm $dir/tmp0_T.sam");
	system("java -Xmx7g -jar $picard SamToFastq I=$dir/tmp_filteredbamT.sam F=$dir/test_tmp_sequence_1.fastq F2=$dir/test_tmp_sequence_2.fastq VALIDATION_STRINGENCY=LENIENT 2>> $dir/HISAT2.log.out");
	system("rm $dir/tmp_filteredbamT.sam");
	
        my $size = `du -s $dir/test_tmp_sequence_2.fastq`;
	my @temp = split /\t/, $size;	
	if ($temp[0] == 0){
		system("$hisat2 -p 8 -x $hisat2_genome -U $dir/test_tmp_sequence_1.fastq -S $dir/hisat2_realign.sam");
	}
	else{
		system("$hisat2 -p 8 -x $hisat2_genome -1 $dir/test_tmp_sequence_1.fastq -2 $dir/test_tmp_sequence_2.fastq -S $dir/hisat2_realign.sam");
	}

	system("rm $dir/test_tmp_sequence_1.fastq");
	system("rm $dir/test_tmp_sequence_2.fastq");
	system("samtools sort -T $dir/tmp $dir/hisat2_realign.sam > $dir/hisat2_realign.bam");
	system("samtools index $dir/hisat2_realign.bam");
	system("rm $dir/hisat2_realign.sam");
	
	system("java -jar $picard MarkDuplicates -I $dir/hisat2_realign.bam -O $dir/md_hisat2_realign.bam -M $dir/md_hisat2_realign.bam_md_metrics.txt -VALIDATION_STRINGENCY STRICT 2>> $dir/HISAT2.log.out");
	system("$gatk SplitNCigarReads -R $fastaDir -I $dir/md_hisat2_realign.bam -O $dir/split_md_hisat2_realign.bam 2>> $dir/HISAT2.log.out");
	system("rm $dir/md_hisat2_realign.bam");
	system("java -jar $picard AddOrReplaceReadGroups -I $dir/split_md_hisat2_realign.bam -O $dir/add_split_md_hisat2_realign.bam -RGID $tumorID -RGLB $tumorLB -RGPL $tumorPL -RGPU $tumorPU -RGSM $tumorSM -VALIDATION_STRINGENCY STRICT 2>> $dir/HISAT2.log.out");
	system("rm $dir/split_md_hisat2_realign.bam");
	system("$gatk BaseRecalibrator -R $fastaDir -I $dir/add_split_md_hisat2_realign.bam -known-sites $dbsnpDir -O $dir/hisat2_realign.table 2>> $dir/HISAT2.log.out");
	system("$gatk ApplyBQSR -R $fastaDir -I $dir/add_split_md_hisat2_realign.bam -bqsr-recal-file $dir/hisat2_realign.table -O $dir/recal_hisat2_realign.bam 2>> $dir/HISAT2.log.out");
	system("rm $dir/add_split_md_hisat2_realign.bam");
	system("rm $dir/hisat2_realign.bam");
	
	my $end_time_epoch = time;
	my $end_time_str = getTime();
	my $first_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str END Extract Reads from $dir Variants Calling\n";
	print "$end_time_str Time Period: ", $first_epoch, "\n";
	print "##################################################################################################################################################\n";
	
}

sub getTime{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	my $time_str = sprintf("%02d-%02d-%02d %02d:%02d:%02d", $year-100,$mon+1, $mday, $hour, $min, $sec);
	return $time_str;
}

sub timeTranslate{
	my ($start,$end) = @_;
	my $epoch = $end - $start;
	use integer;
	$epoch =  sprintf("%02d:%02d:%02d", $epoch/3600, $epoch/60%60, $epoch%60);
	return $epoch;
}
