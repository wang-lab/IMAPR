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
		$input =~ /-ID\t/
		&& $input =~ /-mode\t/
		&& $input =~ /-T\t/
		&& $input =~ /-N\t/		
		&& $input =~ /-O\t/
		
		#system
		&& $input =~ /-thread\t/
		&& $input =~ /-ram\t/
		
		#tools
		&& $input =~ /-gatk\t/
		&& $input =~ /-picard\t/
		&& $input =~ /-hisat2\t/
		&& $input =~ /-samtools\t/
		
		#reference
		&& $input =~ /-R\t/
		&& $input =~ /-gtf\t/
		&& $input =~ /-gene\t/
		&& $input =~ /-dbsnp\t/
		&& $input =~ /-hisat2_reference\t/
		&& $input =~ /-germline\t/
		&& $input =~ /-pon\t/
		)		
		|| $input eq '-h'){
	print "Found missing paramters, Please use -h for help information\n";
	exit;	
}

if ($input eq '-h'){
	print "This is the script for discovering variants based on RNA-seq data. \n";
	print "Usage: perl detect_varaints.pl [options]... -ID [sample ID] -O [output directory] -R [fasta reference]...\n";
	print "###version: 1.0.0\n";
	print "############################################################Parameters############################################################\n";
	print "Required:\n";
	print "##########Input files##########\n";
	print "-ID\n\tsample_name: sample name\n\n";
	print "-mode\n\tinput_format: input files format, RNA/RNA OR RNA/WXS\n\n";
	print "-T\n\ttumor_input: path to aligned bam file for tumor sample\n\n";
	print "-N\n\tnormal_input: path to aligned bam file for normal sample\n\n";
	print "-O\n\tout_prefix: path to output folder\n\n";
	
	print "##########System setting##########\n";
	print "-thread\n\tthreads: number of threads\n\n";
	print "-ram\n\tram: available ram to run IMAPR in gigabytes(GB), 30GB is required to run IMAPR\n\n";
	
	print "##########Tools##########\n";
	print "-gatk\n\tgatk: path to gatk tool package\n\n";
	print "-picard\n\tpicard: path to picard jar file\n\n";
	print "-hisat2\n\thisat2: path to hisat2 package\n\n";
	print "-samtools\n\tsamtools: path to samtools package\n\n";
	
	print "##########Reference##########\n";
	print "-R\n\tfasta_ref: path to genome fasta reference\n\n";
	print "-gtf\n\tgtf_ref: path to gtf reference\n\n";
	print "-gene\n\tgenelist_ref: path to gene list reference\n\n";
	print "-dbsnp\n\tdbsnp_ref: path to dbsnp reference\n\n";
	print "-hisat2_reference\n\thisat_ref: path to hisat2 reference\n\n";
	print "-germline\n\tgermline_ref: path to germline reference\n\n";
	print "-pon\n\tPON_ref: path to PON reference\n\n";
	
	exit;
}

###############################################################################
#prepare the input data from arugument input
my $idInput = '.';
if ($input =~ /-ID\t(\S+)/){
	$idInput = $1;
}

my $modeInput = '.';
if ($input =~ /-mode\t(\S+)/){
	$modeInput = $1;
}

my $tumorInput = '.';
if ($input =~ /-T\t(\S+)/){
	$tumorInput = $1;
}

my $normalInput = '.';
if ($input =~ /-N\t(\S+)/){
	$normalInput = $1;
}

my $fastaReference = '.';
if ($input =~ /-R\t(\S+)/){
	$fastaReference = $1;
}

my $outDir = '.';
if ($input =~ /-O\t(\S+)/){
	$outDir = $1;
}

###############################################################################
my $ram = '.';
if ($input =~ /-ram\t(\S+)/){
	$ram = $1;
}
my $thread = '.';
if ($input =~ /-thread\t(\S+)/){
	$thread = $1;
}
###############################################################################
my $gtfFile = '.';
if ($input =~ /-gtf\t(\S+)/){
	$gtfFile = $1;
}

my $geneList = '.';
if ($input =~ /-gene\t(\S+)/){
	$geneList = $1;
}


my $dbsnpFile = '.';
if ($input =~ /-dbsnp\t(\S+)/){
	$dbsnpFile = $1;
}

my $hisat2Reference = '.';
if ($input =~ /-hisat2_reference\t(\S+)/){
	$hisat2Reference = $1;
}

my $germlineResource = '.';
if ($input =~ /-germline\t(\S+)/){
	$germlineResource = $1;
}

my $ponResource = '.';
if ($input =~ /-pon\t(\S+)/){
	$ponResource = $1;
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

my $samtools = '.';
if ($input =~ /-samtools\t(\S+)/){
	$samtools = $1;
}
#################################################################################################################
print "##################################################################################################################################################\n";
print "##########################################################                              ##########################################################\n";
print "##########################################################  Running detect_variants.pl  ##########################################################\n";
print "##########################################################                              ##########################################################\n";
print "##################################################################################################################################################\n";
#############################################Check RAM requirement###############################################
if ($ram < 30){
	print "Error: Doesn't meet basic requirement for the memory, please try to run IMAPR with additional RAM. \n";
	exit;
}

###########################################Generate output directory#############################################
unless(-d $outDir){
	system("mkdir $outDir");
}

##############################################Preprocess bam file################################################
my %gtf = readGTF($gtfFile);
if($modeInput eq 'RNA/DNA'){
	
	#preprocess for tumor input
	my $inputResource = "Tumor";
	
	#check bam file index
	my $tumorBamIndex = $tumorInput;	
	$tumorBamIndex =~ s/\.bam/\.bai/;
	my $tumorBamIndex2 = $tumorInput . '.bai';	
	unless(-e $tumorBamIndex || -e $tumorBamIndex2){		
		system("$samtools index $tumorInput");
	}
	
	my $normalBamIndex = $normalInput;	
	$normalBamIndex =~ s/\.bam/\.bai/;
	my $normalBamIndex2 = $normalInput . '.bai';
	unless(-e $normalBamIndex || -e $normalBamIndex2){
		system("$samtools index $normalInput");
	}
	
	
	my $tumor_ID_string = `$samtools view -H $tumorInput | grep SM: | head -1`;
	unless ($tumor_ID_string =~ /ID\:/){
		print "Error: Can't find sample ID in your input $tumorInput, please use gatk AddOrReplaceReadGroups to add sample ID. \n";
		exit;
	}
	
	RNAPreProcessParallel($outDir, $tumorInput, $inputResource);
	my $tumorInputMutect = "$outDir/reAligned_$inputResource.bam";
	my $normalInputMutect = "$normalInput";
	
	#1st variant detection	
	my $variationIDInput = $idInput . "_first";
	variantsCalling($tumorInputMutect,$normalInputMutect,$variationIDInput,$outDir,$fastaReference);
	
	#extract Reads from tumor bam and realign with Hisat2
	
	$tumor_ID_string = `$samtools view -H $tumorInputMutect | grep SM: | head -1`;
	my $tumorSM = "NA";
	my $tumorID = $1 if $tumor_ID_string =~ /ID\:(\S+)\s/;
	my $tumorPL = "NA";
	my $tumorLB = "NA";
	my $tumorPU = "NA";
	$tumorSM = $1 if $tumor_ID_string =~ /SM\:(\S+)\s/;
	$tumorPL = $1 if $tumor_ID_string =~ /PL\:(\S+)\s/;
	$tumorLB = $1 if $tumor_ID_string =~ /LB\:(\S+)\s/;
	$tumorPU = $1 if $tumor_ID_string =~ /PU\:(\S+)\s/;
	
	extractReadFromBam($outDir, "$outDir/$variationIDInput\_Variants.bed", $tumorInput, $tumorSM, $tumorID, $tumorPL, $tumorLB, $tumorPU);
	
	#2nd variant detection	
	$variationIDInput = $idInput . "_final";
	my $tumorInputMutectReAlign = "$outDir/reAligned_hisat2_$inputResource.bam";
	variantsCalling($tumorInputMutectReAlign,$normalInputMutect,$variationIDInput,$outDir,$fastaReference);
	system("rm -r $outDir/*bedSplit*");
}
elsif($modeInput eq 'RNA/RNA'){
	
	#preprocess for tumor input
	my $inputResource = "Tumor";
	
	#check bam file index
	my $tumorBamIndex = $tumorInput;	
	$tumorBamIndex =~ s/\.bam/\.bai/;
	my $tumorBamIndex2 = $tumorInput . '.bai';	
	unless(-e $tumorBamIndex || -e $tumorBamIndex2){
		system("$samtools index $tumorInput");
	}
	
	my $normalBamIndex = $normalInput;	
	$normalBamIndex =~ s/\.bam/\.bai/;
	my $normalBamIndex2 = $normalInput . '.bai';	
	unless(-e $normalBamIndex || -e $normalBamIndex2){
		system("$samtools index $normalInput");
	}
	
	my $tumor_ID_string = `$samtools view -H $tumorInput | grep SM: | head -1`;
	unless ($tumor_ID_string =~ /ID\:/){
		print "Error: Can't find sample ID in your input $tumorInput, please use gatk AddOrReplaceReadGroups to add sample ID. \n";
		exit;
	}
	
	RNAPreProcessParallel($outDir, $tumorInput, $inputResource);
	my $tumorInputMutect = "$outDir/reAligned_$inputResource.bam";
	
	#preprocess for normal input
	$inputResource = "Normal";
	my $normal_ID_string = `$samtools view -H $normalInput | grep SM: | head -1`;
	unless ($normal_ID_string =~ /ID\:/){
		print "Error: Can't find sample ID in your input $normalInput, please use gatk AddOrReplaceReadGroups to add sample ID. \n";
		exit;
	}
	RNAPreProcessParallel($outDir, $normalInput, $inputResource);
	my $normalInputMutect = "$outDir/reAligned_$inputResource.bam";
		
	#1st variant detection
	my $variationIDInput = $idInput . "_first";
	variantsCalling($tumorInputMutect,$normalInputMutect,$variationIDInput,$outDir,$fastaReference);
	
	#extract Reads from tumor bam and realign with Hisat2	
	$tumor_ID_string = `$samtools view -H $tumorInputMutect | grep SM: | head -1`;
	my $tumorSM = "NA";
	my $tumorID = $1 if $tumor_ID_string =~ /ID\:(\S+)\s/;
	my $tumorPL = "NA";
	my $tumorLB = "NA";
	my $tumorPU = "NA";
	$tumorSM = $1 if $tumor_ID_string =~ /SM\:(\S+)\s/;
	$tumorPL = $1 if $tumor_ID_string =~ /PL\:(\S+)\s/;
	$tumorLB = $1 if $tumor_ID_string =~ /LB\:(\S+)\s/;
	$tumorPU = $1 if $tumor_ID_string =~ /PU\:(\S+)\s/;
	extractReadFromBam($outDir, "$outDir/$variationIDInput\_Variants.bed", $tumorInput, $tumorSM, $tumorID, $tumorPL, $tumorLB, $tumorPU);
	
	#2nd variant detection	
	$variationIDInput = $idInput . "_final";
	my $tumorInputMutectReAlign = "$outDir/reAligned_hisat2_$inputResource.bam";
	variantsCalling($tumorInputMutectReAlign,$normalInputMutect,$variationIDInput,$outDir,$fastaReference);
	system("rm -r $outDir/*bedSplit*");
}
else{
	print "Error: Can't recognize input_format. The input format should be RNA/RNA or RNA/DNA.\n";
	exit;
}

###################################################Subroutines###################################################
sub readGTF{
	my ($filename) = @_;
	my %gtf;
	print "##################################################################################################################################################\n";	
	my $start_time_epoch = time; #epoch time 	
	my $start_time_str = getTime();
	print "$start_time_str START Loading Reference File for $filename\n";
	
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
			print "Error: Can't find the gene information in $filename, please check your gtf annotation file\n";
			exit;
		} 
	}
	else{
		print "Error: Can't find the GTF annotation file in $filename\n";
		exit;
	}
	
	my $end_time_str = getTime();
	my $end_time_epoch = time;
	my $pre_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	
	print "$end_time_str END Loading Reference File for $filename\n";
	print "$end_time_str Time Period: ", $pre_epoch, "\n";
	print "##################################################################################################################################################\n";	
	return  %gtf;
}

sub RNAPreProcessParallel{
	my ($outDir,$bamFile,$format) = @_;
	
	print "##################################################################################################################################################\n";	
	system("rm $outDir/GATK.log.out") if -e "$outDir/GATK.log.out";
	system("rm -r $outDir/tmp") if -d "$outDir/tmp";
	system("rm -r $outDir/reAligned.bam") if -d "$outDir/reAligned.bam";
	system("rm -r $outDir/reAligned.bam.bai") if -d "$outDir/reAligned.bam.bai";
	my $start_time_epoch = time; #epoch time 	
	my $start_time_str = getTime();
	
	print "$start_time_str START Pre-process for $bamFile\n";
	system("mkdir $outDir/tmp") unless -d "$outDir/tmp";
	system("java -jar $picard SplitSamByNumberOfReads -I $bamFile -O $outDir/tmp -N_READS 10000000 -VALIDATION_STRINGENCY LENIENT 2>> $outDir/GATK.log.out");
	
	my $inter_time_epoch = time;
	my $inter_time_str = getTime();
	print "$inter_time_str FINISH Splitting bam file for $bamFile\n";
	
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
	
	$inter_time_epoch = time;
	$inter_time_str = getTime();
	print "$inter_time_str FINISH Removing Duplicates Reads for $bamFile\n";
	
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
	
	$inter_time_epoch = time;
	$inter_time_str = getTime();
	print "$inter_time_str FINISH Splitting Reads based on Introns for $bamFile\n";
	
	my $BRCommand = "$gatk BaseRecalibrator -R $fastaReference ";
	foreach my $i(1..$fileCount){
		$i = sprintf("%04d", $i);
		$BRCommand = $BRCommand . "-I $outDir/tmp/split_md_shard_$i.bam ";		
	}
	$BRCommand = $BRCommand . "-known-sites $dbsnpFile -O $outDir/tmp/Aligned.out.table 1>> $outDir/GATK.log.out 2>> $outDir/GATK.log.out";	
	system($BRCommand);	
	
	$inter_time_epoch = time;
	$inter_time_str = getTime();
	print "$inter_time_str FINISH Reads Recalibration for $bamFile\n";
	
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
	
	$inter_time_epoch = time;
	$inter_time_str = getTime();
	print "$inter_time_str FINISH BQSR for $bamFile\n";
	
	my $mergeCommand = "java -jar $picard MergeSamFiles ";
	foreach my $i(1..$fileCount){
		$i = sprintf("%04d", $i);
		$mergeCommand = $mergeCommand . "-I $outDir/tmp/re_split_md_shard_$i.bam ";		
	}
	$mergeCommand = $mergeCommand . "-O $outDir/reAligned_$format.bam 2>> $outDir/GATK.log.out";

	system($mergeCommand);	
	$inter_time_epoch = time;
	$inter_time_str = getTime();
	print "$inter_time_str FINISH Merging Pre-processed bam file for $bamFile\n";
	
	system("rm -r $outDir/tmp");
	system("$samtools index $outDir/reAligned_$format.bam");
	
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
	
	my $variantFile = "$outDir/$ID\_Variants.vcf";
	my $variantBed = "$outDir/$ID\_Variants.bed";
	system("rm -r $variantFile") if -e "$variantFile";	
	system("rm -r $variantBed") if -e "$variantBed";	
	
	#test if bed file exist
	
	my $lineCount = `wc -l $geneList`;
	$lineCount = $1 if $lineCount =~ /(\d+)\s/;
	my $lineBreaker = ceil($lineCount/$thread);
		
	system("rm -r $outDir/bedSplit*.list") if -e "$outDir/bedSplit*.list";
	system("split -l $lineBreaker $geneList $outDir/bedSplit");
	my @bedList = `ls $outDir/bedSplit*`;
	

	my $tumor_ID_string = `$samtools view -H $tumorInput | grep SM: | head -1`;
	unless ($tumor_ID_string =~ /ID\:/){
		print "Error: Can't find sample ID in your input $tumorInput, please use gatk AddOrReplaceReadGroups to add sample ID. \n";
		exit;
	}
	my $tumor_ID = $1 if $tumor_ID_string =~ /SM\:(\S+)\s/;
	
	
	my $normal_ID_string = `$samtools view -H $normalInput | grep SM: | head -1`;
	unless ($normal_ID_string =~ /ID\:/){
		print "Error: Can't find sample ID in your input $normalInput, please use gatk AddOrReplaceReadGroups to add sample ID. \n";
		exit;
	}
	my $normal_ID = $1 if $normal_ID_string =~ /SM\:(\S+)\s/;
		
	my $command = '';	
	
	foreach my $bedList (@bedList){
		$bedList =~ s/\s+$//;
		next if $bedList =~ /vcf/;
		next if $bedList =~ /list/;
		system("mv $bedList $bedList.list") unless -e "$bedList.list";		
		my $bed = $1 if $bedList =~ /$outDir\/(\S+)/;
		#my $bedStatus = "failed";
		
		system("rm $outDir/$ID\_$bed.log.out") if -e "$outDir/$ID\_$bed.log.out";
		system("rm $outDir/$ID\_$bed.vcf") if -e "$outDir/$ID\_$bed.vcf";
		system("rm $outDir/$ID\_$bed.vcf.idx") if -e "$outDir/$ID\_$bed.vcf.idx";
		system("rm $outDir/$ID\_$bed.vcf.stats") if -e "$outDir/$ID\_$bed.vcf.stats";
		
		$command = $command . "$gatk Mutect2 -R $fastaReference -I $tumorInput -I $normalInput -tumor $tumor_ID -normal $normal_ID -O $outDir/$ID\_$bed.vcf -L $bedList.list --panel-of-normals $ponResource --germline-resource $germlineResource --read-validation-stringency LENIENT --disable-read-filter AllowAllReadsReadFilter 1>> $outDir/$ID\_$bed.log.out 2>> $outDir/$ID\_$bed.log.out &\n";
	}
	$command = $command . "wait\n";	
	system($command);
	
	foreach my $bedList (@bedList){
		$bedList =~ s/\s+$//;
		system("rm $bedList") if -e $bedList;
		system("rm $bedList.list") if -e "$bedList.list";
	}		
	
	my %vcf;
	my @vcf = `ls $outDir/$ID*.vcf`;
	foreach my $vcf (@vcf){
		$vcf =~ s/\s+$//;
		#print "$vcf\n";
		open(VCF, "$vcf") or die "Cannot open file for reading: $!\n";
		while (<VCF>) {
			next if $_ =~ /^\#/;	
			$_ =~ s/\s+$//;
			my @line = split /\t/, $_;			
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
		close(VCF);
	}
	
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
	print "$start_time_str START Extract Reads from $bamInput 1st Variants Calling\n";
	
	my $tmp_dir = "$dir/extract_tmp";
	my $hisat2_log = "$dir/HISAT2.log.out";
	my $hisat2_stats = "$dir/HISAT2.stats";	
	
	system("rm -r $tmp_dir") if -d "$tmp_dir";
	system("rm $hisat2_log") if -e "$hisat2_log";
	system("rm $hisat2_stats") if -e "$hisat2_stats";	
	system("mkdir $tmp_dir") unless -d "$tmp_dir";		
	
	open(my $IDs_all, ">", "$tmp_dir/IDs_all.txt") or die "Cannot open IDs_all.txt: $!";
	open(my $samtools_view, "-|", "$samtools view -L $variantsFile -@ $thread $bamInput") or die "Cannot run $samtools view: $!";
	while (my $line = <$samtools_view>) {
		my @fields = split("\t", $line);
		print $IDs_all "$fields[0]\n";
	}
	close($samtools_view);
	close($IDs_all);

	system("java -Xmx7g -jar $picard FilterSamReads -I $bamInput -O $tmp_dir/tmp_bam.bam -READ_LIST_FILE $tmp_dir/IDs_all.txt -FILTER includeReadList -WRITE_READS_FILES false -VALIDATION_STRINGENCY LENIENT 2>> $hisat2_log");

	open(my $tmp_filteredbamT, ">", "$tmp_dir/tmp_filteredbamT.sam") or die "Cannot open tmp_filteredbamT.sam: $!";
	open(my $samtools_view_header, "-|", "$samtools view -H $tmp_dir/tmp_bam.bam") or die "Cannot run $samtools view: $!";
	while (my $line = <$samtools_view_header>) {
		print $tmp_filteredbamT $line unless $line =~ /^\@PG/;
	}
	close($samtools_view_header);

	open(my $samtools_view_body, "-|", "$samtools view $tmp_dir/tmp_bam.bam") or die "Cannot run $samtools view: $!";
	while (my $line = <$samtools_view_body>) {
		my @fields = split("\t", $line);
		print $tmp_filteredbamT $line if $fields[1] < 2040;
	}
	close($samtools_view_body);
	close($tmp_filteredbamT);

	system("java -Xmx7g -jar $picard SamToFastq -I $tmp_dir/tmp_filteredbamT.sam -F $tmp_dir/test_tmp_sequence_1.fastq -F2 $tmp_dir/test_tmp_sequence_2.fastq -VALIDATION_STRINGENCY LENIENT 2>> $hisat2_log");
	
	my $size = `du -s $tmp_dir/test_tmp_sequence_2.fastq`;
	my @temp = split /\t/, $size;
	
	if ($temp[0] == 0){
		system("$hisat2 -p 8 -x $hisat2Reference -U $tmp_dir/test_tmp_sequence_1.fastq -S $tmp_dir/hisat2_realign.sam 2>> $hisat2_stats");
	}
	else{
		system("$hisat2 -p 8 -x $hisat2Reference -1 $tmp_dir/test_tmp_sequence_1.fastq -2 $tmp_dir/test_tmp_sequence_2.fastq -S $tmp_dir/hisat2_realign.sam 2>> $hisat2_stats");
	}
	
	my $aligmentOnePer = 0;
	my $aligmentZeroPer = 0;
	my $aligmentMultiPer = 0;
	open(STATS, $hisat2_stats) or die "Cannot open hisat2 alignment stats file $hisat2_stats for reading: $!\n";
	while (<STATS>) {
		
	}
	close(STATS);
	
	system("$samtools sort -T $tmp_dir $tmp_dir/hisat2_realign.sam > $tmp_dir/hisat2_realign.bam 2>/dev/null");
	system("$samtools index $tmp_dir/hisat2_realign.bam");
	
	system("java -jar $picard MarkDuplicates -I $tmp_dir/hisat2_realign.bam -O $tmp_dir/md_hisat2_realign.bam -M $tmp_dir/md_hisat2_realign.bam_md_metrics.txt -VALIDATION_STRINGENCY STRICT 2>> $hisat2_log");
	system("$gatk SplitNCigarReads -R $fastaReference -I $tmp_dir/md_hisat2_realign.bam -O $tmp_dir/split_md_hisat2_realign.bam 2>> $hisat2_log");	
	system("java -jar $picard AddOrReplaceReadGroups -I $tmp_dir/split_md_hisat2_realign.bam -O $tmp_dir/add_split_md_hisat2_realign.bam -RGID $tumorID -RGLB $tumorLB -RGPL $tumorPL -RGPU $tumorPU -RGSM $tumorSM -VALIDATION_STRINGENCY STRICT 2>> $hisat2_log");
	system("$gatk BaseRecalibrator -R $fastaReference -I $tmp_dir/add_split_md_hisat2_realign.bam -known-sites $dbsnpFile -O $tmp_dir/hisat2_realign.table 1>> $hisat2_log 2>> $hisat2_log");
	system("$gatk ApplyBQSR -R $fastaReference -I $tmp_dir/add_split_md_hisat2_realign.bam -bqsr-recal-file $tmp_dir/hisat2_realign.table -O $dir/reAligned_hisat2_Tumor.bam 2>> $hisat2_log");
	
	system("rm -r $tmp_dir") if -d $tmp_dir;
	my $end_time_epoch = time;
	my $end_time_str = getTime();
	my $first_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
	print "$end_time_str END Extract Reads from $bamInput 1st Variants Calling\n";
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
