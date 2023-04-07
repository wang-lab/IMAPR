use strict;
use warnings;
use lib '/home/gtang/data1';
#use XW::Seq;
#use XW::File;
use Data::Dumper;
use List::Util qw/sum/;
use List::Util qw/max/;
use List::Util qw/min/;
#use Statistics::Test::WilcoxonRankSum;
#use Capture::Tiny qw/capture/;

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
my $pon = "/data1/gtang/genome/gdc/MuTect2.PON.5210.vcf";
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



my $tag;

my @cancer = ("GBM");
foreach my $cancer(@cancer){
	print "$cancer\n";
	my @list = `ls -d ./$cancer/TCGA-*/`;
	foreach my $list(@list){
		$list =~ s/\s+$//;
		#print "$list\n";
		my @tcga = split /\//, $list;
		
		my %wgs;
		my %wxs;
=commnet		
		my $wxs = "/data3/gtang/server/$cancer/mutation/mutect_vcf/$tcga[2]/$tcga[2].maf";
		next unless -e $wxs;
		open(IN, "$wxs") or die "Cannot open $wxs for reading: $!\n";	
		while(<IN>){
			next if $_ =~ /^\#/;
			$_=~ s/\s+$//;
			my @line = split /\t/, $_;
			next unless $line[6] eq 'PASS';
			#next unless $line[6] eq 'PASS' or $line[6] eq 'panel_of_normals';
			$wxs{$line[0]}{$line[1]} = $_;
		}
		close(IN);
=cut		
				
		my $vcf = "./$cancer/$tcga[2]/final_Variants.vcf";
		my $vcf_raw = "./$cancer/$tcga[2]/first_Variants.vcf";
		
		next unless -e $vcf;
		next unless -e $vcf_raw;
		
		

		
		my %raw;
		open(IN, "$vcf_raw") or die "Cannot open $vcf_raw for reading: $!\n";
		while(<IN>){
			next if $_ =~ /^\#/;
			$_=~ s/\s+$//;
			my @line = split /\t/, $_;
			
			$raw{$line[0]}{$line[1]} = 1;
			#print "$line[0]\t$line[1]\n";
			
		}
		close(IN);
		
		
		my %vcf;
		open(IN, "$vcf") or die "Cannot open $vcf for reading: $!\n";
		while(<IN>){
			next if $_ =~ /^\#/;
			$_=~ s/\s+$//;
			my @line = split /\t/, $_;
			
			$vcf{$line[0]}{$line[1]} = $_;
			#print "$line[0]\t$line[1]\n";
			
		}
		close(IN);
		my %tmp;
		#print "$tcga[2]\n";
		my $tag = $tcga[2];
		my %pass_vcf;		
		
		#update bed file
		my $bed_file = "./$cancer/$tcga[2]/final_Variants.bed";
		my $update_bed_file = "./$cancer/$tcga[2]/update_final_Variants.bed";
		update_bed($bed_file,$update_bed_file,1,0,'mpileup');
		
		
		my $bam = `ls ./$cancer/$tcga[2]/*/recal_hisat2_realign.bam`;
		$bam =~ s/\s+$//;
		my $mpileupout = "./$cancer/$tcga[2]/mpileup.output";
		
		system ("samtools mpileup -d 1000000 -l $update_bed_file -f $fastq_reference $bam -o $mpileupout 2>/dev/null") unless -e $mpileupout;
		
		$update_bed_file = "./$cancer/$tcga[2]/update_fastq_final_Variants.bed";
		update_bed($bed_file,$update_bed_file,4,4,'fasta');
		my $variant_fasta = "./$cancer/$tcga[2]/variant_fasta.output";
		
		system ("samtools faidx -r $update_bed_file $fastq_reference -o $variant_fasta") unless -e $variant_fasta;
		
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
				my $pileup_out;
				#check if there is any reads aligned in this region by samtools.
				if(exists $mpileup{$chr}{$pos}){
					my ($seq, $pos, $ref_base, $coverage, $pileup, $quals) = split /\t/, $mpileup{$chr}{$pos};
					my @temp =  split /\t/, $mpileup{$chr}{$pos};
					$pileup_out = $pileup;
					#check if there is any bases was detected with alteration in tumor sample by samtools.
					
				
				
						
					
					
					if ($pileup =~ /[ACGTNacgtn]/){
						my $i = 0;
						
							my $temp = length($pileup);
							print "$temp\t";
							$temp = length($quals);
							print "$temp\t";
							print "$mpileup{$chr}{$pos}\n";
						while ($i < length($pileup)) {
							my $char = substr($pileup, $i, 1);

							if ($char eq '^') {
								$i += 2; # Skip the next character (mapping quality)
								next;
							} elsif ($char eq '$') {
								$i++; # Just skip the current character (read end)
								next;
							}
							elsif ($char eq '+' or $char eq '-'){
								$i += 3;
								next;
							}
							
							if ($char eq "."){
								push @mapping_qualities_ref, ord(substr($quals, 0, 1)) - 33;
								push @ref_postive, $char;
								$quals = substr($quals, 1); # Remove the first character from $quals
							}
							elsif($char eq ",") {
								push @mapping_qualities_ref, ord(substr($quals, 0, 1)) - 33;
								push @ref_negative, $char;
								$quals = substr($quals, 1); # Remove the first character from $quals
								
							}
							elsif ($char =~ /[ACGTN]/) {
								push @mapping_qualities_alt, ord(substr($quals, 0, 1)) - 33;
								push @alt_postive, $char;
								$quals = substr($quals, 1); # Remove the first character from $quals
							}
							elsif($char=~ /[acgtn]/){
								push @mapping_qualities_alt, ord(substr($quals, 0, 1)) - 33;
								push @alt_negative, $char;
								$quals = substr($quals, 1); # Remove the first character from $quals
							}
							$i++;
						}
						my $alt_length = scalar(@mapping_qualities_alt);
						$mpileup_flag = 1 if $alt_length < 2;						
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
				
				#my $sequencing_quality_prob = sequencing_quality($chr,$pos,$bam);
				#print "$vcf{$chr}{$pos}\t$sequencing_quality_prob\n";
				
				#rb_filter
				my $rb_flag = 0;
				foreach my $pos1 (keys %{$rb{$chr}}){
					if ($pos > $pos1 and $pos < $rb{$chr}{$pos1}){
						$rb_flag = 1;
					}
				}
				next if $rb_flag == 1;
				
				
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
				#$tumor_flag = 1 if $tumor_af < 0.1;
				#if ($tumor_af > 0.1){
				#	$tumor_flag = 0;
				#}
				#elsif ($tumor_frac > 3){
				#	$tumor_flag = 0;
				#}
				
				#my $mpileupout = mpileup($chr,$pos,$bam);
				#print "$chr\t$pos\t$tumor_depth\t$mpileupout\n";
				
				
				
				my $bcf_coverage;
				#my $rna = `samtools depth -r $chr:$pos-$pos $bam`;
				#my $mpile = `samtools mpileup -r $chr:$pos-$pos -f /data1/gtang/genome/gdc/GRCh38.d1.vd1.fa $bam`;
				
				
				
				$tumor_flag = 0 if $tumor_af >= 0.1 and $tumor_frac >= 3;
				next if $tumor_flag ==1;
				
				#coverage_check
				my $coverage_flag = 0;
				my $fake_flag = 0;
				my $high_flag = 0;
				if($tumor_depth > 10 and $tumor_af > 0.95){
					$high_flag = 1;
				}
				next if $high_flag == 1;
				$vcf_flag = $vcf_flag . "High_AF_in_variants;" if $high_flag == 1;
				$vcf_flag = $vcf_flag . "Variants_in_clipping;" if $coverage_flag == 1 or $fake_flag == 1;
				
				
				
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
				
				my $Mono_flag = 0;
				
				my $l1 = scalar @alt_postive;
				my $l2 = scalar @alt_negative;
				#print "$pileup_out\t$l1\t$l2\n";
				if(length($vcf_flag) ==0){					
					my @seq = split//, $fasta{$chr}{$pos};
					#print "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$sub[1]\n";
					if(length($line[3]) == 1 and length($line[4]) == 1){
						$seq[4] = $line[4];
						my $sub2 = join( '', @seq);
						if($line[4] eq $seq[3] or $line[4] eq $seq[5]){
							my $change = "$line[4]$line[4]$line[4]$line[4]";
							if($sub2 =~ $change){
								#$Mono_flag = 1 if ($l1+1)/($l2+1) > 2 or ($l1+1)/($l2+1) < 0.5;
								$Mono_flag = 1;
								#print "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$sub[1]\n";
								
							}
						}
					}				
				}
				next if $Mono_flag ==1;
				$vcf_flag = $vcf_flag . "Tandem_Repeats;" if $Mono_flag == 1;
		
				
				
				
=comment			
				if(length($vcf_flag) ==0){
					my $bam = `ls ./$cancer/$tcga[2]/*/recal_hisat2_realign.bam`;
					#system("samtools index $bam");
					$bam =~ s/\s+$//;
					#print "$bam\n";
					my $start = $line[1] - 2;
					my $end = $line[1] + 2;
					my $rna = `samtools depth -r $line[0]:$start-$end $bam`;
					my @rna = split /\n/, $rna;
					my @coverage = (0,0,0,0,0);
					#print "$line[0]\t$line[1]\n";
					#print "$rna\n";
					if(length($rna) == 0){
						@coverage = (0,0,0,0,0);
					}
					else{
						foreach my $rna1(@rna){
							my @read = split /\t/, $rna1;
							if($read[1] eq $start){
								$coverage[0] = $read[2];
							}
							elsif($read[1] eq ($start+1)){
								$coverage[1] = $read[2];
							}
							elsif($read[1] eq $end){
								$coverage[4] = $read[2];
							}
							elsif($read[1] eq ($end-1)){
								$coverage[3] = $read[2];
							}
							elsif($read[1] eq $line[1]){
								$coverage[2] = $read[2];
							}
						}
					}
					my $coverage = join( ',', @coverage );
					
					#print "@coverage\n";
					if ($coverage[2] == 0){
						$coverage_flag = 1;
					}
					elsif(max(@coverage)/$coverage[2] > 1.5){
						$coverage_flag = 1;
					}
					elsif(($coverage[2]+1)/(min(@coverage)+1) > 1.5){
						$coverage_flag = 1;
					}
					my $ratio = 0;
					if ($coverage[2] != 0){
						$ratio = $tumor_depth/$coverage[2];
						$fake_flag = 1 if $ratio > 1.5;
					}

				}
=cut				
				
				#fake reads
				
				
				
				#$vcf_flag = "";
				#$vcf_flag = $vcf_flag . "panel_of_normal;" if $pon_flag == 1;
				$vcf_flag = "PASS" if length($vcf_flag) == 0;
				$pass_vcf{$chr}{$pos} = $vcf{$chr}{$pos};
				#$pass_vcf{$chr}{$pos} = $vcf{$chr}{$pos} if $vcf_flag eq "PASS";				
				
			}
		}
		print "$tcga[2]\t";
		compare(\%pass_vcf,\%wxs,\%wgs,\%tmp);
		
	
	}
}



sub compare{
	my ($rna,$wxs,$wgs,$tmp) = @_;
	my %rna = %$rna;
	my %wxs = %$wxs;
	my %wgs = %$wgs;
	my %tmp = %$tmp;
	my $total_count = 0;
	my $share_count = 0;
	my $wxs_count = 0;
	my $wgs_count = 0;
	
	
	
	foreach my $chr(keys %rna){
		foreach my $pos(keys %{$rna{$chr}}){
			#print "$tag\t";
			$total_count++;
			if( (exists $wxs{$chr}{$pos} || exists $wxs{$chr}{$pos+1} || exists $wxs{$chr}{$pos+2} ||  exists $wxs{$chr}{$pos+3} || exists $wxs{$chr}{$pos+4} || exists $wxs{$chr}{$pos+5} || exists $wxs{$chr}{$pos-1} || exists $wxs{$chr}{$pos-2} || exists $wxs{$chr}{$pos-3} || exists $wxs{$chr}{$pos-4} || exists $wxs{$chr}{$pos-5} ) 
			&& (exists $wgs{$chr}{$pos} || exists $wgs{$chr}{$pos-1} || exists $wgs{$chr}{$pos-2} || exists $wgs{$chr}{$pos-3} || exists $wgs{$chr}{$pos-4} || exists $wgs{$chr}{$pos-5} || exists $wgs{$chr}{$pos+1} || exists $wgs{$chr}{$pos+2} || exists $wgs{$chr}{$pos+3} || exists $wgs{$chr}{$pos+4} || exists $wgs{$chr}{$pos+5} ) ){
				$share_count++;
				#print OUT "$tmp{$chr}{$pos}\t1\n";
				#print "$tmp{$chr}{$pos}\t";
				#print "wxs\twgs\t$rna{$chr}{$pos}\n";
			}
			elsif(exists $wxs{$chr}{$pos} || exists $wxs{$chr}{$pos+1} || exists $wxs{$chr}{$pos+2} ||  exists $wxs{$chr}{$pos+3} || exists $wxs{$chr}{$pos+4} || exists $wxs{$chr}{$pos+5} || exists $wxs{$chr}{$pos-1} || exists $wxs{$chr}{$pos-2} || exists $wxs{$chr}{$pos-3} || exists $wxs{$chr}{$pos-4} || exists $wxs{$chr}{$pos-5}){
				$wxs_count++;
				#print OUT "$tmp{$chr}{$pos}\t1\n";
				#print "$tmp{$chr}{$pos}\t";
				#print "wxs\tno\t$rna{$chr}{$pos}\n";
			}
			elsif(exists $wgs{$chr}{$pos} || exists $wgs{$chr}{$pos-1} || exists $wgs{$chr}{$pos-2} || exists $wgs{$chr}{$pos-3} || exists $wgs{$chr}{$pos-4} || exists $wgs{$chr}{$pos-5} || exists $wgs{$chr}{$pos+1} || exists $wgs{$chr}{$pos+2} || exists $wgs{$chr}{$pos+3} || exists $wgs{$chr}{$pos+4} || exists $wgs{$chr}{$pos+5}){
				$wgs_count++;
				#print OUT "$tmp{$chr}{$pos}\t1\n";
				#print "$tmp{$chr}{$pos}\t";
				#print "no\twgs\t$rna{$chr}{$pos}\n";
			}
			else{
				#print OUT "$tmp{$chr}{$pos}\t0\n";
				#print "$tmp{$chr}{$pos}\t";
				#print "no\tno\t$rna{$chr}{$pos}\n";
			}
		}
	}
	print "$total_count\t$share_count\t$wxs_count\t$wgs_count\n";
	
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
	}
	close(IN);
	close(OUT);
}

sub mpileup{
	my ($ref_seq,$position,$bam_file) = @_;
	my ($mpileup_output, $stderr) = capture {
		system ( "samtools mpileup -r ${ref_seq}:${position}-${position} -f /data1/gtang/genome/gdc/GRCh38.d1.vd1.fa $bam_file" );
	};
	$mpileup_output =~ s/\s+$//;
	return $mpileup_output;	
}

sub sequencing_quality{
	my ($ref_seq,$position,$bam_file) = @_;
	
	my $flag = '';
	my ($mpileup_output, $stderr) = capture {
		system ( "samtools mpileup -r ${ref_seq}:${position}-${position} -f /data1/gtang/genome/gdc/GRCh38.d1.vd1.fa $bam_file" );
	};	
	
	$mpileup_output =~ s/\s+$//;
	my @mapping_qualities_ref;
	my @mapping_qualities_alt;

	
	
	my ($seq, $pos, $ref_base, $coverage, $pileup, $quals) = split /\t/, $mpileup_output;
	#print "$mpileup_output\n";
	if ($mpileup_output eq ''){
		$flag = $flag . "Low_reads;";
	}
	elsif (length($quals) <= 1){
		$flag = $flag . "Low_reads;";
	} 
	elsif ($pileup !~ /[ACGTNacgtn]/){
		$flag = $flag . "Non_alteration_reads;";
	}
	else{
		my $i = 0;
		while ($i < length($pileup)) {
			my $char = substr($pileup, $i, 1);

			if ($char eq '^') {
				$i += 2; # Skip the next character (mapping quality)
				next;
			} elsif ($char eq '$') {
				$i++; # Just skip the current character (read end)
				next;
			}

			if ($char eq "." || $char eq ",") {
				push @mapping_qualities_ref, ord(substr($quals, 0, 1)) - 33;
				$quals = substr($quals, 1); # Remove the first character from $quals
			} elsif ($char =~ /[ACGTNacgtn]/) {
				push @mapping_qualities_alt, ord(substr($quals, 0, 1)) - 33;
				$quals = substr($quals, 1); # Remove the first character from $quals
			}

			$i++;
		}
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
		$wilcox_test->load_data(\@mapping_qualities_ref, \@mapping_qualities_alt);
		my $prob = $wilcox_test->probability();
		#$flag = $flag . $prob;
	}
	
	return $flag;
}

