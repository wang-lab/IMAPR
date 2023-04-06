use strict;
use warnings;
use Data::Dumper;
use POSIX;



#tools
my $gatk = "/home/mcho43/tools/gatk-4.1.8.1/gatk-4.1.8.1/gatk";
my $picard = "/home/mcho43/tools/gatk-4.1.8.1/gatk-4.1.8.1/picard.jar";
my $hisat2 = "/home/mcho43/tools/hisat2-2.1.0/hisat2-2.1.0/hisat2";
my $gdcClient = "/home/mcho43/tools/gdc-client";
my $gdcToken = "/home/mcho43/tools/gdc-user-token.2023-03-13T21_44_46.673Z.txt";
my $star = "STAR";

#indices
my $genomeDir = "/projects/pharma_xwang317/mcho43/genome";
my $STARGenomeDir = "/projects/pharma_xwang317/mcho43/genome";
my $fastaDir = "/projects/pharma_xwang317/mcho43/genome/GRCh38.d1.vd1.fa";
my $gtfDir = "/projects/pharma_xwang317/mcho43/genome/gencode.v36.annotation.gtf";
my $dbsnpDir = "/projects/pharma_xwang317/mcho43/genome/common_all_20180418_update.vcf.gz";
my $gtfBedDir = "/projects/pharma_xwang317/mcho43/genome/gg.list";
my $hisat2_genome = "/projects/pharma_xwang317/mcho43/genome/hisat2_GRCh38";

#param
my $thread = "12";
my $ram = "96";
my $sampleID;
my $samplePL;
my $sampleLB;
my $samplePU;
my $patID = $ARGV[0];


my $gtfFlag = 1;
my %gtf;
if ($gtfFlag == 1){
	open(GTF, "$gtfDir") or die "Cannot open $gtfDir for reading: $!\n";
	while (<GTF>) {
		next if $_ =~ /^\#/;
		$_ =~ s/\s+$//;
		my @line = split /\t/, $_;
		next unless $line[2] eq 'exon';
		$gtf{$line[0]}{$line[3]} = 1;
		$gtf{$line[0]}{$line[4]} = 1;
	}
	close(GTF);
}
my $outHead = "/projects/pharma_xwang317/mcho43/TCGA";
my $m_file = "$outHead/$cancer_type/gdc_m.txt";
open(IN, "$m_file") or die "Cannot open $m_file for reading: $!\n";
while (<IN>) {
	#next if $. < 2;
	next unless $. eq $patID;
	$_ =~ s/\s+$//;
	my @line = split /\t/, $_;
	my $command = "mkdir $outHead/$cancer_type/$line[0]";
	system($command) unless -d "$outHead/$cancer_type/$line[0]";	
	
	$command = "$gdcClient download -t $gdcToken $line[3] -d $outHead/$cancer_type/$line[0] -n $thread";
	system($command);
	
	my $tumorDir = "$outHead/$cancer_type/$line[0]/$line[3]";
	my $normalDir = "$outHead/$cancer_type/$line[0]/$line[1]";
	my $tcgaDir = "$outHead/$cancer_type/$line[0]";
	my $tumorBam = "$tumorDir/$line[4]";
	my $normalBam = "$normalDir/$line[2]";
	
	my $tumor_ID_string = `samtools view -H $tumorBam | grep SM: | head -1`;
	my $tumorSM = $1 if $tumor_ID_string =~ /SM\:(\S+)\s/;
	my $tumorID = $1 if $tumor_ID_string =~ /ID\:(\S+)\s/;
	my $tumorPL = $1 if $tumor_ID_string =~ /PL\:(\S+)\s/;
	my $tumorLB = $1 if $tumor_ID_string =~ /LB\:(\S+)\s/;
	my $tumorPU = $1 if $tumor_ID_string =~ /PU\:(\S+)\s/;
	
	print "$tumorSM\n$tumorID\n$tumorPL\n$tumorLB\n$tumorPU\n";
	
	#RNAPreProcess($tumorDir,$tumorDir,$line[4]);
	RNAPreProcessParallel($tumorDir,$tumorDir,$line[4]);
	
	$command = "$gdcClient download -t $gdcToken $line[1] -d $outHead/$cancer_type/$line[0] -n $thread";
	system($command);

	my $tumorPreBam = "$tumorDir/reAligned.bam";
	my $variantsCallingID = "first";
        
        my $normalBamIndex = $normalBam;	
	$normalBamIndex =~ s/\.bam/\.bai/;
	my $normalBamIndex2 = $normalBam . '.bai';	
	unless(-e $normalBamIndex or $normalBamIndex2){
		system("samtools index $normalBam");
	}
	variantsCalling($tumorPreBam,$normalBam,$variantsCallingID,$tcgaDir,$genomeDir);
	
	my $extractBam = "$tumorDir/$line[4]";
	extractReadFromBam($tumorDir, "$tcgaDir/first_Variants.bed", $extractBam, $tumorSM, $tumorID, $tumorPL, $tumorLB, $tumorPU);
	$variantsCallingID = "final";
	$tumorBam = "$tumorDir/recal_hisat2_realign.bam";
	variantsCalling($tumorBam,$normalBam,$variantsCallingID,$tcgaDir,$genomeDir);
	system("rm $tumorPreBam");
	system("rm $normalBam");
	system("rm $tumorDir/$line[4]");
	
}
close(IN);

sub RNAPreProcessParallel{
	my ($inputFile,$outDir,$bamFile) = @_;
	my $start_time_epoch = time; #epoch time 
	
	my $start_time_str = getTime();
	print "##################################################################################################################################################\n";
	print "$start_time_str START RNA-seq Aligment for $inputFile\n";
	 
	system("rm $outDir/STAR.log.out");
	system("rm $outDir/GATK.log.out");
	#system("$star --readFilesIn $inputFile --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimOutType Junctions SeparateSAMold --chimSegmentMin 15 --genomeDir $genomeDir --limitSjdbInsertNsj 1200000 --outFileNamePrefix $outDir --outFilterIntronMotifs None --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.66 --outFilterType BySJout --outSAMattributes NH HI AS nM NM ch --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --runThreadN $thread --twopassMode Basic --outSAMattrRGline ID:$sampleID PL:$samplePL LB:$sampleLB PU:$samplePU SM:tumor --outSAMmapqUnique 60 >> $outDir/STAR.log.out");
	my $end_time_epoch = time;
	my $end_time_str = getTime();
	print "$end_time_str END RNA-seq Aligment for $inputFile\n";
	my $align_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str Time Period: ", $align_epoch, "\n";
	
	$start_time_str = getTime();
	$start_time_epoch = time;
	print "$start_time_str START Pre-process for $inputFile\n";
	system("mkdir $outDir/tmp") unless -d "$outDir/tmp";
	system("java -jar $picard SplitSamByNumberOfReads -I $outDir/$bamFile -O $outDir/tmp -N_READS 10000000 -VALIDATION_STRINGENCY LENIENT 2>> $outDir/GATK.log.out");
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
		#print "$markDupCommand\n";
		system($markDupCommand);
		
	}
	
	$upLimit = floor($ram/10);
	#print "$upLimit\n";
	$iteration = ceil($fileCount/$upLimit);
	foreach my $i(1..$iteration){
		my $splitCommand = "";
		foreach my $j(1..$upLimit){
			my $temp = ($i-1)*$upLimit+$j;
			last if $temp > $fileCount;
			$temp = sprintf("%04d", $temp);
			$splitCommand = $splitCommand . "$gatk SplitNCigarReads -R $fastaDir -I $outDir/tmp/md_shard_$temp.bam -O $outDir/tmp/split_md_shard_$temp.bam 2>> $outDir/GATK.log.out &\n";
		}
		$splitCommand = $splitCommand . "wait\n";
		#print "$splitCommand\n";
		system($splitCommand);
		
	}
	my $BRCommand = "$gatk BaseRecalibrator -R $fastaDir ";
	foreach my $i(1..$fileCount){
		$i = sprintf("%04d", $i);
		$BRCommand = $BRCommand . "-I $outDir/tmp/split_md_shard_$i.bam ";		
	}
	$BRCommand = $BRCommand . "-known-sites $dbsnpDir -O $outDir/tmp/Aligned.out.table 2>> $outDir/GATK.log.out";
	#print "$BRCommand\n";
	system($BRCommand);	
	
	foreach my $i(1..$iteration){
		my $BQSRCommand = "";
		foreach my $j(1..$upLimit){
			my $temp = ($i-1)*$upLimit+$j;
			last if $temp > $fileCount;
			$temp = sprintf("%04d", $temp);
			$BQSRCommand = $BQSRCommand . "$gatk ApplyBQSR -R $fastaDir -I $outDir/tmp/split_md_shard_$temp.bam -bqsr-recal-file $outDir/tmp/Aligned.out.table -O $outDir/tmp/re_split_md_shard_$temp.bam 2>> $outDir/GATK.log.out &\n";
		}
		$BQSRCommand = $BQSRCommand . "wait\n";
		#print "$BQSRCommand\n";
		system($BQSRCommand);
		
	}
	
	my $mergeCommand = "java -jar $picard MergeSamFiles ";
	foreach my $i(1..$fileCount){
		$i = sprintf("%04d", $i);
		$mergeCommand = $mergeCommand . "-I $outDir/tmp/re_split_md_shard_$i.bam ";		
	}
	$mergeCommand = $mergeCommand . "-O $outDir/reAligned.bam 2>> $outDir/GATK.log.out";
	#print "$mergeCommand\n";
	system($mergeCommand);
	system("rm -r $outDir/tmp");
	system("samtools index $outDir/reAligned.bam");
	
	$end_time_str = getTime();
	$end_time_epoch = time;
	my $pre_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str END Pre-process for $inputFile\n";
	print "$end_time_str Time Period: ", $pre_epoch, "\n";
	print "##################################################################################################################################################\n";
}

sub RNAPreProcess{
	my ($inputFile,$outDir,$bamFile) = @_;
	my $start_time_epoch = time; #epoch time 
	
	my $start_time_str = getTime();
	print "##################################################################################################################################################\n";
	print "$start_time_str START RNA-seq Aligment for $inputFile\n";
	
	system('rm $outDir/STAR.log.out');
	system('rm $outDir/GATK.log.out');
	
	#system("$star --readFilesIn $inputFile --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimOutType Junctions SeparateSAMold --chimSegmentMin 15 --genomeDir $genomeDir --limitSjdbInsertNsj 1200000 --outFileNamePrefix $outDir --outFilterIntronMotifs None --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.66 --outFilterType BySJout --outSAMattributes NH HI AS nM NM ch --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --runThreadN $thread --twopassMode Basic --outSAMattrRGline ID:$sampleID PL:$samplePL LB:$sampleLB PU:$samplePU SM:tumor --outSAMmapqUnique 60 >> $outDir/STAR.log.out");
	my $end_time_epoch = time;
	my $end_time_str = getTime();
	print "$end_time_str END RNA-seq Aligment for $inputFile\n";
	my $align_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str Time Period: ", $align_epoch, "\n";
	
	$start_time_str = getTime();
	$start_time_epoch = time;
	print "$start_time_str START Pre-process for $inputFile\n";
	system("java -jar $picard MarkDuplicates -I $outDir/$bamFile -O $outDir/md_Aligned.out.bam -M $outDir/md_Aligned.out.bam_md_metrics.txt --REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT 2>> $outDir/GATK.log.out");
	system("rm $outDir/$bamFile");
	system("$gatk SplitNCigarReads -R $fastaDir -I $outDir/md_Aligned.out.bam -O $outDir/split_md_Aligned.out.bam 2>> $outDir/GATK.log.out");
	system("rm $outDir/md_Aligned.out.bam");
	system("$gatk BaseRecalibrator -R $fastaDir -I $outDir/split_md_Aligned.out.bam -known-sites $dbsnpDir -O $outDir/Aligned.out.table 2>> $outDir/GATK.log.out");
	system("$gatk ApplyBQSR -R $fastaDir -I $outDir/split_md_Aligned.out.bam -bqsr-recal-file $outDir/Aligned.out.table -O $outDir/recal_Aligned.out.bam 2>> $outDir/GATK.log.out");
	system("rm $outDir/split_md_Aligned.out.bam");
	
	$end_time_str = getTime();
	$end_time_epoch = time;
	my $pre_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str END Pre-process for $inputFile\n";
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
