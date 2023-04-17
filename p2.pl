use strict;
use warnings;
use lib '/home/gtang/data1';
#use XW::Seq;
#use XW::File;
use Data::Dumper;
use List::Util qw/sum/;
use List::Util qw/max/;
use List::Util qw/min/;
use Statistics::Test::WilcoxonRankSum;
use Capture::Tiny qw/capture/;
use Math::CDF qw(pnorm);
use Statistics::RankCorrelation;
use Statistics::Distributions;
use Math::Complex;
use Statistics::Distributions qw(chisqrprob);

#########################This is the script for filtering variants based on RNA-seq data#########################
########################The input file should be an output directory generated from p1.pl########################

my $start_time_epoch = time; #epoch time 
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $start_time_str = sprintf("%02d-%02d-%02d_%02d-%02d-%02d", $year-100,$mon+1, $mday, $hour, $min, $sec);

my $input = join "\t", @ARGV;
$input =~ s/\s+$//;
unless ((
		#input_files
		$input =~ /-iputDir\t/		
		)		
		|| $input eq '-h'){
	print "Found missing paramters, Please use -h for help information\n";
	exit;	
}

if ($input eq '-h'){
	print "This is the script for filtering variants based on RNA-seq data. The input file should be an output directory generated from p1.pl.\n";
	print "Usage: perl [options]... -mode single -d ./data_files -o ./out_files -g ./GRCH38 -f pair-end\n";
	print "The result will be stored in output directory, name with Final_*";
	print "###version: 1.0.0\n";
	print "#########################################################################################################################################\n";
	print "Required:\n";	
	print "-iputDir\n\t./data_files: The directory for input file\n\n";	
	exit;
}

###############################################################################
#prepare the input data from arugument input
my $inputDir = '.';
if ($input =~ /-iputDir\t(\S+)/){
	$inputDir = $1;
}


###############################################################################
my $fastq_reference = "/data1/gtang/genome/gdc/GRCh38.d1.vd1.fa";

my %igg;
my $ref = "/data1/gtang/genome/gdc/IGG.txt";
open(IN, "$ref") or die "Cannot open $ref for reading: $!\n";
while(<IN>){
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$igg{$line[0]}{$line[3]} = $line[4];
}
close(IN);

my %hla;
$ref = "/data1/gtang/genome/gdc/HLA.txt";
open(IN, "$ref") or die "Cannot open $ref for reading: $!\n";
while(<IN>){
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$hla{$line[0]}{$line[3]} = $line[4];
}
close(IN);

my %rb;
$ref = "/data1/gtang/genome/gdc/ribosomal_protein.txt";
open(IN, "$ref") or die "Cannot open $ref for reading: $!\n";
while(<IN>){
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$rb{$line[0]}{$line[3]} = $line[4];
}
close(IN);

my %pseudo;
$ref = "/data1/gtang/genome/gdc/pseudo_gene.txt";
open(IN, "$ref") or die "Cannot open $ref for reading: $!\n";
while(<IN>){
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$pseudo{$line[0]}{$line[3]} = $line[4];
}
close(IN);

my %pon;
my %ponDetail;
my $pon = "/data1/gtang/ref/MuTect2.PON.5210.vcf";
open(IN, "$pon") or die "Cannot open $pon for reading: $!\n";
while (<IN>) {
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	#next if length($line[3]) > 5;
	my @target = split /\,/, $line[4];
	foreach my $target(@target){
		#next if length($target) > 5;
		if(length($line[3]) ==1 and length($target) == 1){
			$ponDetail{$line[0]}{$line[1]}{"$line[3]_$target"} = 1;
		}
	}
	$pon{$line[0]}{$line[1]} = 1;
}
close(IN);

my %edits;
my $edits = "/data1/gtang/genome/TABLE1_hg38.txt";
open(IN, "$edits") or die "Cannot open $pon for reading: $!\n";
while (<IN>) {
	next if $. < 2;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$edits{$line[0]}{$line[1]} = 1;
	$edits{$line[0]}{$line[1]+1} = 1;
	$edits{$line[0]}{$line[1]-1} = 1;
}
close(IN);

$edits = "/data1/gtang/genome/gdc/Darned_38.bed";
open(IN, "$edits") or die "Cannot open $pon for reading: $!\n";
while (<IN>) {
	#next if $. < 2;
	$_=~ s/\s+$//;
	my @line = split /[:-]/, $_;
	$edits{$line[0]}{$line[1]} = 1;
	$edits{$line[0]}{$line[1]+1} = 1;
	$edits{$line[0]}{$line[1]-1} = 1;
}
close(IN);

$edits = "/data1/gtang/genome/gdc/Radar_38.bed";
open(IN, "$edits") or die "Cannot open $pon for reading: $!\n";
while (<IN>) {
	#next if $. < 2;
	$_=~ s/\s+$//;
	my @line = split /[:-]/, $_;
	$edits{$line[0]}{$line[1]} = 1;
	$edits{$line[0]}{$line[1]+1} = 1;
	$edits{$line[0]}{$line[1]-1} = 1;
}
close(IN);
	
my $vcf = "./$cancer/$tcga[2]/final_Variants.vcf";
my $vcf_raw = "./$cancer/$tcga[2]/first_Variants.vcf";
		
my %raw;
open(IN, "$vcf_raw") or die "Cannot open $vcf_raw for reading: $!\n";
while(<IN>){
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;			
	$raw{$line[0]}{$line[1]} = 1;
}
close(IN);
		
		
my %vcf;
open(IN, "$vcf") or die "Cannot open $vcf for reading: $!\n";
while(<IN>){
	next if $_ =~ /^\#/;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;			
	$vcf{$line[0]}{$line[1]} = $_;	
}
close(IN);

my %tmp;
my $tag = $tcga[2];
my %pass_vcf;		
		
#update bed file
my $bed_file = "./$cancer/$tcga[2]/final_Variants.bed";
my $update_bed_file = "./$cancer/$tcga[2]/update_final_Variants.bed";
update_bed($bed_file,$update_bed_file,1,0,'mpileup');
		
my $bam = `ls ./$cancer/$tcga[2]/*/recal_hisat2_realign.bam`;
$bam =~ s/\s+$//;
my $mpileupout = "./$cancer/$tcga[2]/mpileup2.output";		
system ("samtools mpileup -d 1000000 -l $update_bed_file -f $fastq_reference $bam -o $mpileupout 2>/dev/null") unless -e $mpileupout;
		
$update_bed_file = "./$cancer/$tcga[2]/update_fastq_final_Variants.bed";
update_bed($bed_file,$update_bed_file,8,8,'fasta');
my $variant_fasta = "./$cancer/$tcga[2]/variant_fasta.output";		
system ("samtools faidx -r $update_bed_file $fastq_reference -o $variant_fasta") unless -e $variant_fasta;
		
$update_bed_file = "./$cancer/$tcga[2]/update_depth_final_Variants.bed";
update_bed($bed_file,$update_bed_file,3,2,'depth');
my $depthOut = "./$cancer/$tcga[2]/depthOut.output";
system ("samtools depth -b $update_bed_file $bam -o $depthOut") unless -e $depthOut;
		
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
		
my @list = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY");
foreach my $chr(@list){
	foreach my $pos(sort {$a <=> $b} keys %{$vcf{$chr}}){
				my $vcf_flag = "";
				
				my @line = split /\t/, $vcf{$chr}{$pos};
				#duo_alignment
				my $duo_flag = 1;
				if (exists $raw{$chr}{$pos}){				
					$duo_flag = 0;
				}				
				next if $duo_flag == 1;
				
				#cluster_filter
				my $cluster_flag = 0;			
				if ($line[8] =~ /PGT/){
					$cluster_flag = 1;			
				}
				next if $cluster_flag == 1;
				
				#edits_filter
				my $edits_flag = 0;			
				if (exists $edits{$chr}{$pos}){
					$edits_flag = 1;			
				}
				next if $edits_flag == 1;
								
				#igg_filter
				my $igg_flag = 0;
				foreach my $pos1 (keys %{$igg{$chr}}){
					if ($pos > $pos1 and $pos < $igg{$chr}{$pos1}){
						$igg_flag = 1;
					}
				}
				next if $igg_flag == 1;
				
				#hla_filter
				my $hla_flag = 0;
				foreach my $pos1 (keys %{$hla{$chr}}){
					if ($pos > $pos1 and $pos < $hla{$chr}{$pos1}){
						$hla_flag = 1;
					}
				}
				next if $hla_flag == 1;
				
				#pon_filter
				my $pon_flag = 0;			
				#if (exists $pon{$chr}{$pos} or exists $pon{$chr}{$pos+1} or exists $pon{$chr}{$pos-1}){
				#	$pon_flag = 1;			
				#}
				if(length($line[3]) == 1 and length($line[4]) == 1){
					if (exists $pon{$chr}{$pos} and exists $ponDetail{$chr}{$pos}{"$line[3]_$line[3]"}){
						$pon_flag = 1;
					}
				}
				else{
					if (exists $pon{$chr}{$pos} or exists $pon{$chr}{$pos+1} or exists $pon{$chr}{$pos-1}){
						$pon_flag = 1;			
					}
				}
				next if $pon_flag == 1;
				
				#length_filter
				my $length_flag = 0;			
				if (length($line[3]) > 5 or length($line[4]) > 5){
					$length_flag = 1;			
				}
				#next if $length_flag == 1;
				
				#normal_cutoff
				my $normal_flag = 0;
				my @normal = split /:/, $line[10];
				my @normal_af = split /,/, $normal[1];
				my @normal_per = split /,/, $normal[2];
				my $normal_af = max(@normal_per);
				shift @normal_af;
				my $normal_frac = max(@normal_af);
				#next if $normal_af > 0.1;		##high AF
				#next if $normal_frac > 1;		##high fraction
				if ($normal_af > 0.1 or $normal_frac > 1){
					$normal_flag = 1;
				}
				next if $normal_flag == 1;
				
				#PON
				if ($line[7] =~ /PON/){
					$pon_flag = 1;			
				}
				next if $pon_flag == 1;
				
				#tandem_repeats
				my $tandem_flag = 0;
				if ($line[7] =~ /STR/ or $line[7] =~ /RU/){
					$tandem_flag = 1;
				}
				next if $tandem_flag == 1;
				
				#multiallelic
				my $multiallelic_flag = 0;			
				if ($line[4] =~ /\,/){
					$multiallelic_flag = 1;			
				}
				next if $multiallelic_flag == 1;
				
				#mpileup
				my $mpileup_flag = 0;
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
				if(exists $mpileup{$chr}{$pos}){
					my ($seq, $pos, $ref_base, $coverage, $pileup, $quals) = split /\t/, $mpileup{$chr}{$pos};
					$pileup_out = $pileup;
					#check if there is any bases was detected with alteration in tumor sample by samtools.
					my $temp_char = '';
					if ($pileup =~ /[ACGTNacgtn]/){
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
						#my $alt_length = scalar(@mapping_qualities_alt);
						
						$mpileup_flag = 1 if ($forward_alt+$reverse_alt) <= 2;	
						$mpileup_flag = 1 if ($forward_alt+$reverse_alt)/($forward_all+$reverse_all) < 0.1;	
					}
					else{
						$mpileup_flag = 1;
					}
				}
				else{
					$mpileup_flag = 1;
				}
				next if $mpileup_flag == 1;
				
				#sequencing quality
			
				
				
				#rb_filter
				my $rb_flag = 0;
				foreach my $pos1 (keys %{$rb{$chr}}){
					if ($pos > $pos1 and $pos < $rb{$chr}{$pos1}){
						$rb_flag = 1;
					}
				}
				#next if $rb_flag == 1;
				
				
				#pseudo_filter
				my $pseudo_flag = 0;
				foreach my $pos1 (keys %{$pseudo{$chr}}){
					if ($pos > $pos1 and $pos < $pseudo{$chr}{$pos1}){
						$pseudo_flag = 1;
					}
				}
				next if $pseudo_flag == 1;
				
				
				#
				my $MPOS = 100;
				$MPOS = $1 if $line[7] =~ /;MPOS=(\d+);/ ;
				my $pos_flag = 0;
				if($MPOS < 6){
					$pos_flag = 1;
				}
				next if $pos_flag ==1;
				
				my $TLOD = 100;
				$TLOD = $1 if $line[7] =~ /;TLOD=(\S+)/;
				$TLOD = 100 if $TLOD =~ /\,/; 
				my $TLOD_flag = 0;
				if($TLOD < 5.6){
					$TLOD_flag = 1;					
				}
				next if $TLOD_flag == 1;
				$tmp{$chr}{$pos} = $TLOD;
				
				
				my $MMQ = 100;
				$MMQ = $1 if $line[7] =~ /;MMQ=(\S+);MPOS/;
				my $MMQ_flag = 0;
				my @MMQ = split /,/, $MMQ;
				foreach my $mmq (@MMQ){
					if($mmq < 50){
						$MMQ_flag = 1;
					}
				}
				next if $MMQ_flag ==1;
				
				
				#tumor_cutoff
				my $tumor_flag = 1;
				my @tumor = split /:/, $line[9];			
				my @tumor_af = split /,/, $tumor[1];
				my @tumor_per = split /,/, $tumor[2];			
				
				my $tumor_depth =sum(@tumor_af);
				shift @tumor_af;			
				my $tumor_frac = max(@tumor_af);
				my $tumor_af = max(@tumor_per);
				
				$tumor_flag = 0 if $tumor_af >= 0.1 and $tumor_frac >= 3;
				next if $tumor_flag ==1;
				
				#high mutation frequency	
				my $high_flag = 0;
				if($tumor_depth > 10 and $tumor_af > 0.95){
					$high_flag = 1;
				}
				next if $high_flag == 1;
				$vcf_flag = $vcf_flag . "High_AF_in_variants;" if $high_flag == 1;
				
				
				
				
			
							
				my $Mono_flag = 0;				
				my @seq = split//, $fasta{$chr}{$pos};
				my $left_consec;
				my $right_consec;				
			
				#Check single nucleotide changes
				if(length($line[3]) == 1 and length($line[4]) == 1){
					my $mutation_position = 8;
					my ($left_count, $right_count) = count_consecutive_altered_nucleotides($fasta{$chr}{$pos}, $mutation_position, $line[4]);
					if (($left_count+$right_count) >= 4){
						$Mono_flag = 1;
					}
				}				
				next if $Mono_flag ==1;				
				$vcf_flag = $vcf_flag . "Tandem_Repeats;" if $Mono_flag == 1;
				
				
				
			
				my $sequencing_quality_flag = 0;
				my $logFC_sequencing_quality = 0;
				if ($mpileup_flag == 0 and scalar(@mapping_qualities_ref) > 2 and scalar(@mapping_qualities_alt) > 2){
					my $p = mann_whitney_u_test(\@mapping_qualities_ref,\@mapping_qualities_alt);
					$p = $p/2;
					my $med = median(@mapping_qualities_alt);
					my $med2 = median(@mapping_qualities_ref);
					$logFC_sequencing_quality = log(median(@mapping_qualities_alt)/median(@mapping_qualities_ref))/log(2);					
					$sequencing_quality_flag = 1 if $p <= 0.05 and $logFC_sequencing_quality < -0.5;					
				}
				next if $sequencing_quality_flag == 1;
				
				my $coverage_flag = 0;
				#extract the coverage distribution around a mutation site
				my @coverage;
				my $i = $line[1]-2;							
				while ($i <= $line[1]+2) {
					if(exists  $depth{$line[0]}{$i}){
						push @coverage, $depth{$line[0]}{$i};
					}
					else{
						push @coverage, 0;
					}
					$i++;
				}			
				
				my $pvalue = chi_square_test_uniform(@coverage);				
				if($pvalue < 0.05){
					if ($coverage[2] == 0){
						$coverage_flag = 1;
					}
					elsif(max(@coverage)/$coverage[2] > 1.5){
							$coverage_flag = 1;
					}
					elsif(($coverage[2]+1)/(min(@coverage)+1) > 1.5){
						$coverage_flag = 1;
					}					
				}				
				next if $coverage_flag == 1;
				
				$vcf_flag = $vcf_flag . "Single_variants;" if $duo_flag == 1;
				$vcf_flag = $vcf_flag . "Igg_gene;" if $igg_flag == 1;
				$vcf_flag = $vcf_flag . "Hla_gene;" if $hla_flag == 1;
				$vcf_flag = $vcf_flag . "Rb_gene;" if $rb_flag == 1;
				$vcf_flag = $vcf_flag . "Pseudo_gene;" if $pseudo_flag == 1;
				$vcf_flag = $vcf_flag . "panel_of_normal;" if $pon_flag == 1;
				$vcf_flag = $vcf_flag . "RNA_edits;" if $edits_flag == 1;
				$vcf_flag = $vcf_flag . "Long_edits;" if $length_flag == 1;
				$vcf_flag = $vcf_flag . "Cluster_events;" if $cluster_flag == 1;
				$vcf_flag = $vcf_flag . "Germline_risk;" if $normal_flag == 1;
				$vcf_flag = $vcf_flag . "Tandem_Repeats;" if $tandem_flag == 1;
				$vcf_flag = $vcf_flag . "multiallelic;" if $multiallelic_flag == 1;
				$vcf_flag = $vcf_flag . "Softclipping;" if $pos_flag == 1;			
				$vcf_flag = $vcf_flag . "Low_tumor_reads;" if $tumor_flag == 1;
				$vcf_flag = $vcf_flag . "Low_TLOD;" if $TLOD_flag == 1;
				$vcf_flag = $vcf_flag . "Low_mapping_quality;" if $MMQ_flag == 1;					
				
				$vcf_flag = "PASS" if length($vcf_flag) == 0;
				$pass_vcf{$chr}{$pos} = $vcf{$chr}{$pos};
				#$pass_vcf{$chr}{$pos} = $vcf{$chr}{$pos} if $vcf_flag eq "PASS";				
				
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

sub update_bed{
	my ($bed_file,$update_bed_file,$left,$right,$type) = @_;	
	open(IN, "$bed_file") or die "Cannot open $bed_file for reading: $!\n";
	open(OUT, ">$update_bed_file") or die "Cannot open $update_bed_file for reading: $!\n";
	while(<IN>){		
		$_=~ s/\s+$//;
		my @line = split /\t/, $_;
		$line[1] = $line[1] - $left;
		$line[2] = $line[2] + $right;
		print OUT "$line[0]\t$line[1]\t$line[2]\n" if $type eq 'mpileup';
		print OUT "$line[0]:$line[1]-$line[2]\n" if $type eq 'fasta';
		print OUT "$line[0]\t$line[1]\t$line[2]\n" if $type eq 'depth';
	}
	close(IN);
	close(OUT);
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

sub calculate_p_value {
    my ($u, $n1, $n2) = @_;
    my $mean_u = $n1 * $n2 / 2;
    my $var_u = $n1 * $n2 * ($n1 + $n2 + 1) / 12;
    my $z = ($u - $mean_u) / sqrt($var_u);
    my $p_value = 2 * (1 - pnorm(abs($z))); # Two-tailed test
    return $p_value;
}

sub count_consecutive_altered_nucleotides {
    my ($sequence, $position, $altered_base) = @_;
    my $new_sequence = substr($sequence, 0, $position) . $altered_base . substr($sequence, $position + 1);
    
    my $left_seq = substr($new_sequence, 0, $position);
    my $right_seq = substr($new_sequence, $position + 1);
    
    my ($left_count) = $left_seq =~ /($altered_base+)$/; 
    my ($right_count) = $right_seq =~ /^($altered_base+)/;
    
    return (length($left_count || ''), length($right_count || ''));
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