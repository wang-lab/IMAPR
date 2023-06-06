use strict;
use warnings;
use List::Util qw/sum/;
use List::Util qw/max/;
use List::Util qw/min/;
use Statistics::Test::WilcoxonRankSum;
use Statistics::Distributions qw(chisqrprob);

#########################This is the script for filtering variants based on RNA-seq data#########################
########################The input file should be an output directory generated from p1.pl########################

my $input = join "\t", @ARGV;
$input =~ s/\s+$//;
unless ((
		#input_files
		$input =~ /-ID\t/
		&& $input =~ /-O\t/		
		&& $input =~ /-gtf\t/
		)		
		|| $input eq '-h'){
	print "Found missing paramters, Please use -h for help information\n";
	exit;	
}

if ($input eq '-h'){
	print "This is the script for filtering variants using machine learning model .\n";
	print "Usage: perl [options]... -ID [sample ID] -O [output directory] -gtf [gtf reference]\n";	
	print "###version: 1.0.0\n";
	print "#########################################################################################################################################\n";
	print "Required:\n";	
	print "-ID\n\tsample_name: sample name\n\n";
	print "-O\n\tout_prefix: path to output folder\n\n";
	print "-gtf\n\tgtf_ref: path to gtf reference\n\n";
	exit;
}

###############################################################################
#prepare the input data from arugument input
my $idInput = '.';
if ($input =~ /-ID\t(\S+)/){
	$idInput = $1;
}

my $outDir = '.';
if ($input =~ /-O\t(\S+)/){
	$outDir = $1;
}

my $gtfFile = '.';
if ($input =~ /-gtf\t(\S+)/){
	$gtfFile = $1;
}

#################################################################################################################
print "##################################################################################################################################################\n";
print "#########################################################                              ###########################################################\n";
print "#########################################################  Running machine_learning.pl ###########################################################\n";
print "#########################################################                              ###########################################################\n";
print "##################################################################################################################################################\n";

###########################################Generate output directory######################################################
unless(-d $outDir){
	system("mkdir $outDir");
}

##########################################################################################################################

my $start_time_epoch = time; #epoch time
my $start_time_str = getTime();
print "##################################################################################################################################################\n";
print "$start_time_str START Loading Reference File for Machine Learning Model\n";

my %gene;
my %exon;
my %cds;
my %UTR;
if(-e $gtfFile){
	open(GTF, "$gtfFile") or die "Cannot open $gtfFile for reading: $!\n";
	while (<GTF>) {
		next if $_ =~ /^\#/;
		$_ =~ s/\s+$//;
		my @line = split /\t/, $_;
		my $gene_type;
		if ($_ =~ /gene_type "([^"]+)"/) {
			$gene_type = $1;			
		}
		if($gene_type eq 'protein_coding'){
			if($line[2] eq 'gene'){
				if (exists $gene{$line[0]}{$line[3]}){
					$gene{$line[0]}{$line[3]} = $line[4] if $line[4] >= $gene{$line[0]}{$line[3]};
				}
				else{
					$gene{$line[0]}{$line[3]} = $line[4];
				}
			}
			if($line[2] eq 'exon'){
				if (exists $exon{$line[0]}{$line[3]}){
					$exon{$line[0]}{$line[3]} = $line[4] if $line[4] >= $exon{$line[0]}{$line[3]};
				}
				else{
					$exon{$line[0]}{$line[3]} = $line[4];
				}			
			}
			if($line[2] eq 'CDS'){
				if (exists $cds{$line[0]}{$line[3]}){
					$cds{$line[0]}{$line[3]} = $line[4] if $line[4] >= $cds{$line[0]}{$line[3]};
				}
				else{
					$cds{$line[0]}{$line[3]} = $line[4];
				}			
			}
			if($line[2] eq 'UTR'){
				if (exists $UTR{$line[0]}{$line[3]}){
					$UTR{$line[0]}{$line[3]} = $line[4] if $line[4] >= $UTR{$line[0]}{$line[3]};
				}
				else{
					$UTR{$line[0]}{$line[3]} = $line[4];
				}	
			}
		}		
	}
	close(GTF);
	unless (keys %gene or keys %exon or keys %cds or keys %UTR) {
		print "Error: Can't find the gene information in $gtfFile, please check your gtf annotation file\n";
		exit;
	} 
}
else{
	print "Error: Can't find the GTF annotation file in $gtfFile\n";
	exit;
}

my $end_time_epoch = time;
my $end_time_str = getTime();
my $first_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
print "$end_time_str END Loading Reference File for Machine Learning Model\n";
print "$end_time_str Time Period: ", $first_epoch, "\n";
print "##################################################################################################################################################\n";
############################
$start_time_epoch = time; #epoch time
$start_time_str = getTime();
print "##################################################################################################################################################\n";
print "$start_time_str START Running Machine Learning Model\n";

my $out_file = "$outDir/$idInput\_mc_inputs.txt";
open(VAL, ">$out_file") or die "Cannot open $out_file for reading: $!\n";

print VAL "\ttumorDepth\talterTumorReads\talterAF\tmapQuality\tmapPosition\tartifactsNormal\tgenoNormal\tpopulationAF\ttumorLog\twilcoxon\tSQ_alt\tSQ_ref\tSQ_fc\t";
print VAL "mpileup_alt\tmpileup_per\tchi_pvalue\tcoverge_1\tcoverge_2\tcoverge_3\tcoverge_4\tcoverge_5\t";
print VAL "VDB\tSGB\tRPB\tMQB\tMQSB\tBQB\tMQ0F\t";
print VAL "C>T\tC>G\tC>A\tT>G\tT>C\tT>A\t";
print VAL "cds\tUTR\tnc_exon\tintron\tintergenic";
print VAL "\n";

my $vcf = "$outDir/$idInput\_filter_Variants.vcf";
my %vcf;
open(IN, "$vcf") or die "Cannot open $vcf for reading: $!\n";
while(<IN>){
	
	next if $_ =~ /^\#/;	
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$vcf{$line[0]}{$line[1]} = $_;
	
}
close(IN);

my $mpileupout = "$outDir/$idInput\_mpileup.output";
my $bcftools_outfiles = "$outDir/$idInput\_bcftools.output";
my $variant_fasta = "$outDir/$idInput\_variant_fasta.output";
my $depthOut = "$outDir/$idInput\_depthOut.output";

my %mpileup;
open(IN, "$mpileupout") or die "Cannot open $mpileupout for reading: $!\n";
while(<IN>){
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$mpileup{$line[0]}{$line[1]} = $_;
}
close(IN);

my %fasta;
open(IN, "$variant_fasta") or die "Cannot open $variant_fasta for reading: $!\n";
my $temp;
while(<IN>){
	$_=~ s/\s+$//;
	if ($. %2 == 1){
	$temp = $_;
	}
	else{
		my @line = split /[>:-]/, $temp;
		my $med = ($line[2]+$line[3])/2;
		$fasta{$line[1]}{$med} = $_;
	}	
}
close(IN);

my %depth;
open(IN, "$depthOut") or die "Cannot open $depthOut for reading: $!\n";
while(<IN>){
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$depth{$line[0]}{$line[1]} = $line[2];
}
close(IN);

my %bcf;
open(IN, "$bcftools_outfiles") or die "Cannot open $bcftools_outfiles for reading: $!\n";
while(<IN>){
	next if $_ =~ /^#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$bcf{$line[0]}{$line[1]} = $_;
}
close(IN);

my %summary;
foreach my $chr(sort {lc $a cmp lc $b} keys %vcf){
	foreach my $pos(sort {$a <=> $b} keys %{$vcf{$chr}}){		

		my @line = split /\t/, $vcf{$chr}{$pos};
		my @tumor = split /:/, $line[9];			
		my @tumor_af = split /,/, $tumor[1];
		my @tumor_per = split /,/, $tumor[2];			
		
		my $tumor_depth =sum(@tumor_af);
		
		my @mapping_qualities_ref;
		my @mapping_qualities_alt;
		my @ref_postive;
		my @ref_negative;
		my @alt_postive;
		my @alt_negative;
		my $forward_all = 0;
		my $reverse_all = 0;
		my $forward_alt = 0;
		my $reverse_alt = 0;
		my $pileup_out;
		#check if there is any reads aligned in this region by samtools.		
		my ($seq, $pos, $ref_base, $coverage, $pileup, $quals) = split /\t/, $mpileup{$chr}{$pos};
		$pileup_out = $pileup;
		#check if there is any bases was detected with alteration in tumor sample by samtools.
		my $temp_char = '';			
		my $i = 0;
		while ($i < length($pileup)) {
			my $char = substr($pileup, $i, 1);
			if ($char eq '^') {
				$i += 2; # Skip the next character (mapping quality)
				$temp_char = $char;
				next;
			} elsif ($char eq '$') {
				$i++; # Just skip the current character (read end)
				$temp_char = $char;
				next;
			}
			elsif ($char eq '+' or $char eq '-'){
				$i++;
				if($temp_char eq '.'){
					$forward_alt++;
				}
				elsif($temp_char eq ','){
					$reverse_alt++;
				}
				$temp_char = $char;
				next;
			}
			elsif ($char =~ /[0-9]/){
				$i = $i + $char + 1;				
				$temp_char = $char;
				next;
			}
			$temp_char = $char;
			if ($char eq "."){
				push @mapping_qualities_ref, ord(substr($quals, 0, 1)) - 33;
				push @ref_postive, $char;
				$quals = substr($quals, 1); # Remove the first character from $quals
				$forward_all++;
			}
			elsif($char eq ",") {
				push @mapping_qualities_ref, ord(substr($quals, 0, 1)) - 33;
				push @ref_negative, $char;
				$quals = substr($quals, 1); # Remove the first character from $quals
				$reverse_all++;
			}
			elsif ($char =~ /[ACGTN]/) {
				push @mapping_qualities_alt, ord(substr($quals, 0, 1)) - 33;
				push @alt_postive, $char;
				$quals = substr($quals, 1); # Remove the first character from $quals
				$forward_alt++;
				$forward_all++;
			}
			elsif($char=~ /[acgtn]/){
				push @mapping_qualities_alt, ord(substr($quals, 0, 1)) - 33;
				push @alt_negative, $char;
				$quals = substr($quals, 1); # Remove the first character from $quals
				$reverse_alt++;
				$reverse_all++;
			}
			$i++;
		}

		my $p = 0;
		my $med = 0;
		my $med2 = 0;
		my $logFC_sequencing_quality = 0;
		if (scalar(@mapping_qualities_ref) > 2 and scalar(@mapping_qualities_alt) > 2){
			$p = mann_whitney_u_test(\@mapping_qualities_ref,\@mapping_qualities_alt);
			$p = $p/2;
			$logFC_sequencing_quality = log(median(@mapping_qualities_alt)/median(@mapping_qualities_ref))/log(2);
		}
		$med = median(@mapping_qualities_alt) if scalar(@mapping_qualities_alt) > 2;
		$med2 = median(@mapping_qualities_ref) if scalar(@mapping_qualities_ref) > 2;
		$summary{$chr}{$pos}{wil_pvalue} = $p;
		$summary{$chr}{$pos}{seq_med} = $med;
		$summary{$chr}{$pos}{seq_med2} = $med2;
		$summary{$chr}{$pos}{seq_fc} = $logFC_sequencing_quality;
		$summary{$chr}{$pos}{mpileup_alt} = $forward_alt+$reverse_alt;
		$summary{$chr}{$pos}{mpileup_per} = ($forward_alt+$reverse_alt)/($forward_all+$reverse_all);
		
		my @coverage;
		$i = $pos-2;
		while ($i <= $pos+2) {
			if(exists  $depth{$chr}{$i}){
				push @coverage, $depth{$chr}{$i};
			}
			else{
				push @coverage, 0;
			}
			$i++;
		}		
		my $pvalue = chi_square_test_uniform(@coverage);	
		$summary{$chr}{$pos}{chi_pvalue} = $pvalue;
		
		my $maximum = max(@coverage);
		
		for my $i (0 .. $#coverage) {
			$coverage[$i] = $coverage[$i]/$maximum;
		}
		
		$summary{$chr}{$pos}{coverge_1} = $coverage[0];
		$summary{$chr}{$pos}{coverge_2} = $coverage[1];
		$summary{$chr}{$pos}{coverge_3} = $coverage[2];
		$summary{$chr}{$pos}{coverge_4} = $coverage[3];
		$summary{$chr}{$pos}{coverge_5} = $coverage[4];
		
		#VDB=0.0259155;SGB=-0.683931;RPB=0.987467;MQB=1;MQSB=1;BQB=0.999896;MQ0F=0
		my @bcf = split /\t/, $bcf{$chr}{$pos};
		$summary{$chr}{$pos}{VDB} = 0;
		$summary{$chr}{$pos}{SGB} = -1;
		$summary{$chr}{$pos}{RPB} = 0;
		$summary{$chr}{$pos}{MQB} = 0;
		$summary{$chr}{$pos}{MQSB} = 0;
		$summary{$chr}{$pos}{BQB} = 0;
		$summary{$chr}{$pos}{MQ0F} = 0;
		
		$summary{$chr}{$pos}{VDB} = $1 if $bcf[7] =~ /VDB=([\d.]+)/;
		$summary{$chr}{$pos}{SGB} = $1 if $bcf[7] =~ /SGB=(-?[\d.]+)/;
		$summary{$chr}{$pos}{RPB} = $1 if $bcf[7] =~ /RPB=([\d.]+)/;			
		$summary{$chr}{$pos}{MQB} = $1 if $bcf[7] =~ /MQB=([\d.]+)/;			
		$summary{$chr}{$pos}{MQSB} = $1 if $bcf[7] =~ /MQSB=([\d.]+)/;			
		$summary{$chr}{$pos}{BQB} = $1 if $bcf[7] =~ /BQB=([\d.]+)/;			
		$summary{$chr}{$pos}{MQ0F} = $1 if $bcf[7] =~ /MQ0F=([\d.]+)/;			
		
	}
}
generateMCInput(\%vcf,\%summary);	
close VAL;

#machine learning filter
my $mc_out_file = "$outDir/$idInput\_mc_outputs.txt";
system("python3 stacking_model.py $out_file $mc_out_file");

my %mc;
open(MC, "$mc_out_file") or die "Cannot open $mc_out_file for reading: $!\n";
while(<MC>){
	next if $_ =~ /^#/;
	$_=~ s/\s+$//;			
	my @line = split /\t/, $_;
	$mc{$line[1]}{$line[2]} = 1 if $line[3] eq '1';
}
close(MC);

my $mc_filter_file = "$outDir/$idInput\_mc_filter_Variants.vcf";
open(IN, "$vcf") or die "Cannot open $vcf for reading: $!\n";
open(OUT, ">$mc_filter_file") or die "Cannot open $mc_filter_file for reading: $!\n";
while(<IN>){
	next if $_ =~ /^#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	print OUT "$_\n" if exists $mc{$line[0]}{$line[1]};
}
close(IN);
close(OUT);

$end_time_epoch = time;
$end_time_str = getTime();
$first_epoch = timeTranslate($start_time_epoch, $end_time_epoch);
print "$end_time_str END Running Machine Learning Model\n";
print "$end_time_str Time Period: ", $first_epoch, "\n";
print "##################################################################################################################################################\n";

sub generateMCInput{
	my ($rna,$summary) = @_;
	my %rna = %$rna;	
	my %summary = %$summary;
	
	foreach my $chr(sort {lc $a cmp lc $b} keys %rna){
		foreach my $pos(sort {$a <=> $b} keys %{$rna{$chr}}){
					
			my @line = split /\t/, $rna{$chr}{$pos};
			my $id = "$idInput\_$line[0]_$line[1]";
			print VAL "$id\t";			
			my @tumor = split /:/, $line[9];			
			my @tumor_af = split /,/, $tumor[1];
			my @tumor_per = split /,/, $tumor[2];
			
			my $tumor_depth =sum(@tumor_af);
			print VAL "$tumor_depth\t$tumor_af[1]\t$tumor_per[0]\t";
			
			my $MMQ = $1 if $line[7] =~ /MMQ\=(\S+)\;MPOS/;
			my @MMQ = split /\,/, $MMQ;
			print VAL "$MMQ[0]\t";
			
			my $MPOS = $1 if $line[7] =~ /MPOS\=(\d+)\;NALOD/;			
			print VAL "$MPOS\t";			
			
			my $NALOD = $1 if $line[7] =~ /NALOD\=(\S+)\;NLOD/;
			if ($NALOD =~ /e/){
				print VAL "0\t";
			}
			else{
				print VAL "$NALOD\t";
			}
						
			my $NLOD = $1 if $line[7] =~ /NLOD\=(\S+)\;POPAF/;
			print VAL "$NLOD\t";
			
			my $POPAF = $1 if $line[7] =~ /POPAF\=(\S+)\;TLOD/;
			print VAL "$POPAF\t";
			
			my $TLOD = $1 if $line[7] =~ /TLOD\=(\S+)/;
			print VAL "$TLOD\t";
						
			print VAL "$summary{$chr}{$pos}{wil_pvalue}\t";
			print VAL "$summary{$chr}{$pos}{seq_med}\t";
			print VAL "$summary{$chr}{$pos}{seq_med2}\t";
			print VAL "$summary{$chr}{$pos}{seq_fc}\t";
			print VAL "$summary{$chr}{$pos}{mpileup_alt}\t";
			print VAL "$summary{$chr}{$pos}{mpileup_per}\t";
			print VAL "$summary{$chr}{$pos}{chi_pvalue}\t";
			print VAL "$summary{$chr}{$pos}{coverge_1}\t";
			print VAL "$summary{$chr}{$pos}{coverge_2}\t";
			print VAL "$summary{$chr}{$pos}{coverge_3}\t";
			print VAL "$summary{$chr}{$pos}{coverge_4}\t";
			print VAL "$summary{$chr}{$pos}{coverge_5}\t";
			print VAL "$summary{$chr}{$pos}{VDB}\t";
			print VAL "$summary{$chr}{$pos}{SGB}\t";
			print VAL "$summary{$chr}{$pos}{RPB}\t";
			print VAL "$summary{$chr}{$pos}{MQB}\t";
			print VAL "$summary{$chr}{$pos}{MQSB}\t";
			print VAL "$summary{$chr}{$pos}{BQB}\t";
			print VAL "$summary{$chr}{$pos}{MQ0F}\t";			
		
			#nucleotide
			if(($line[3] eq "C" and $line[4] eq "T") or ($line[3] eq "G" and $line[4] eq "A")){
				print VAL "1\t0\t0\t0\t0\t0\t";
			}
			elsif(($line[3] eq "C" and $line[4] eq "G") or ($line[3] eq "G" and $line[4] eq "C")){
				print VAL "0\t1\t0\t0\t0\t0\t";
			}
			elsif(($line[3] eq "C" and $line[4] eq "A") or ($line[3] eq "G" and $line[4] eq "T")){
				print VAL "0\t0\t1\t0\t0\t0\t";
			}
			elsif(($line[3] eq "T" and $line[4] eq "G") or ($line[3] eq "A" and $line[4] eq "C")){
				print VAL "0\t0\t0\t1\t0\t0\t";
			}
			elsif(($line[3] eq "T" and $line[4] eq "C") or ($line[3] eq "A" and $line[4] eq "G")){
				print VAL "0\t0\t0\t0\t1\t0\t";
			}
			elsif(($line[3] eq "T" and $line[4] eq "A") or ($line[3] eq "A" and $line[4] eq "T")){
				print VAL "0\t0\t0\t0\t0\t1\t";
			}
			else{
				print VAL "0\t0\t0\t0\t0\t0\t";
			}
			
			#genomic location
			my $cds_flag = 0;
			foreach my $pos1 (keys %{$cds{$chr}}){
					if ($pos >= $pos1 and $pos <= $cds{$chr}{$pos1}){						
						$cds_flag = 1;
						
						last;
					}
			}
			
			my $utr_flag = 0;
			foreach my $pos1 (keys %{$UTR{$chr}}){
					if ($pos >= $pos1 and $pos <= $UTR{$chr}{$pos1}){
						$utr_flag = 1;
						last;
					}
			}
			
			my $exon_flag = 0;
			foreach my $pos1 (keys %{$exon{$chr}}){
					if ($pos >= $pos1 and $pos <= $exon{$chr}{$pos1}){
						$exon_flag = 1;
						last;
					}
			}
			
			my $gene_flag = 0;
			foreach my $pos1 (keys %{$gene{$chr}}){
					if ($pos >= $pos1 and $pos <= $gene{$chr}{$pos1}){
						$gene_flag = 1;
						last;
					}
			}
			
			if($gene_flag == 1){
				if ($exon_flag == 1){
					if ($cds_flag == 1){
						print VAL "1\t0\t0\t0\t0";
					}
					else{
						if ($utr_flag == 1){
							print VAL "0\t1\t0\t0\t0";
						}
						else{
							print VAL "0\t0\t1\t0\t0";
						}
					}
				}
				else{
					print VAL "0\t0\t0\t1\t0";
				}
			}
			else{
				print VAL "0\t0\t0\t0\t1";
			}			

			print VAL "\n";
			
		}
	}
	
}


sub chi_square_test_uniform {
    my @data = @_;
    my $num_bins = scalar @data;
    my $total = sum @data;
    my $expected_count = $total / $num_bins;
    
    my $chi_square = 0;
    for my $count (@data) {
        my $deviation = $count - $expected_count;
        $chi_square += ($deviation ** 2) / $expected_count;
    }

    my $degrees_of_freedom = $num_bins - 1;
    my $p_value = chisqrprob($degrees_of_freedom, $chi_square);

    return $p_value;
}

sub mann_whitney_u_test {
    my ($array1, $array2) = @_;	
	
	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	$wilcox_test->load_data(\@$array1, \@$array2);
	my $prob = $wilcox_test->probability();
 
	my $pf = sprintf '%f', $prob;
 
	
	my $pstatus = $wilcox_test->probability_status();
	
 
	#$wilcox_test->summary();
	
	
	return $pf;
	
   
}

sub median {
    my @numbers = sort { $a <=> $b } @_;
    my $count = @numbers;

    if ($count % 2) {
        # Odd number of elements, return the middle one
        return $numbers[int($count / 2)];
    } else {
        # Even number of elements, return the average of the two middle ones
        return ($numbers[int($count / 2) - 1] + $numbers[int($count / 2)]) / 2;
    }
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