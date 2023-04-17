use strict;
use warnings;
use lib '/home/gtang/data1';
#use XW::Seq;
#use XW::File;
use Data::Dumper;
use List::Util qw/sum/;
use List::Util qw/max/;
use List::Util qw/min/;

			


open(IN, "/data1/gtang/server/LUAD/signature.txt") or die "Cannot open signature.txt for reading: $!\n";
my %sig;
my %sig_total;
my %pat_total;
while(<IN>){
	next if $. < 2;
	$_=~ s/\s+$//;
	my @line = split /\t/, $_;
	$pat_total{$line[2]} = 1;
	if(exists $sig{$line[1]}{$line[4]}){
		$sig{$line[1]}{$line[4]} = $sig{$line[1]}{$line[4]} + 1;		
	}
	else{
		$sig{$line[1]}{$line[4]} = 1;
	}
	
	if(exists $sig_total{$line[1]}){
		$sig_total{$line[1]} = $sig_total{$line[1]} + 1;		
	}
	else{
		$sig_total{$line[1]} = 1;
	}
}
close(IN);






my @list = `ls -d ./TCGA-*/`;


my $file = "./gdc_m.txt";
open(MA, "$file") or die "Cannot open $file for reading: $!\n";
my $tag;
while (<MA>) {
	next if $. < 2;
	$_ =~ s/\s+$//;
	my @tcga = split /\t/, $_;

	#next unless $. >=$start and $. <= $end;
	my $flag = '';
	foreach my $key (@list){
		my @tmp = split/\//, $key;
		$flag = "PASS" if $tmp[1] eq $tcga[0];
	}
	next unless $flag eq "PASS";
	#next unless $tcga[0] eq "TCGA-55-A493";
	
	my $vcf = "./$tcga[0]/anno.hg38_multianno.txt";
	my $out_file = "./$tcga[0]/RF_inputs.txt";
	open(OUT, ">$out_file") or die "Cannot open $out_file for reading: $!\n";
	print OUT "\ttumorDepth\talterTumorReads\talterAF\tmapQuality\tmapPosition\tartifactsNormal\tgenoNormal\tpopulationAF\ttumorLog\ttypeExonic\ttypeIntronic\ttypeUTR\ttypeStream\ttypeIntergenic\ttypeSplicing\ttypeElse\tmutCounts\tmutRatio\tmutGeneRatio\tclass\n";

	my %wgs;
	my %wxs;	
	
	
	
	
	
	my %vcf;
	open(IN, "$vcf") or die "Cannot open $vcf for reading: $!\n";
	while(<IN>){
		next if $_ =~ /^\#/;
		next if $. < 2;
		$_=~ s/\s+$//;
		my @line = split /\t/, $_;
		
		$vcf{"chr$line[0]"}{$line[1]} = $_;
		#print "$line[0]\t$line[1]\n";
		
	}
	close(IN);
	my %tmp;
	#print "$tcga[0]\t";
	$tag = $tcga[0];
	my %pass_vcf;	
	foreach my $chr(sort {lc $a cmp lc $b} keys %vcf){
		foreach my $pos(sort {$a <=$b} keys %{$vcf{$chr}}){
			my $vcf_flag = "";
			
			my @line = split /\t/, $vcf{$chr}{$pos};
			my @tumor = split /:/, $line[122];			
			my @tumor_af = split /,/, $tumor[1];
			my @tumor_per = split /,/, $tumor[2];			
			
			my $tumor_depth =sum(@tumor_af);
			
			$vcf_flag = $line[119];
			
			#next if $tumor_depth < 10;
			$pass_vcf{$chr}{$pos} = $vcf{$chr}{$pos} if $vcf_flag eq "PASS";
			#print "$tcga[0]\t$vcf{$chr}{$pos}\n" if $vcf_flag eq "PASS";
			
		}
	}
	#print "$tcga[0]\t$wxs_count\t$wgs_count\t";
	compare(\%pass_vcf,\%wxs,\%wgs,\%tmp);
	
	close(OUT);
	
	system("python3 machinelearning.py $out_file ./$tcga[0]");
}
close(MA);



#compare(\%rna, \%wxs, \%wgs);
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
			my @line = split /\t/, $rna{$chr}{$pos};
			my $id = "$tag\_chr$line[113]_$line[114]";
			print OUT "$id\t";
			my @tumor = split /:/, $line[122];			
			my @tumor_af = split /,/, $tumor[1];
			my @tumor_per = split /,/, $tumor[2];			
			
			my $tumor_depth =sum(@tumor_af);
			print OUT "$tumor_depth\t$tumor_af[1]\t$tumor_per[0]\t";
			
			my $MMQ = $1 if $line[120] =~ /MMQ\=(\S+)\;MPOS/;
			my @MMQ = split /\,/, $MMQ;
			print OUT "$MMQ[0]\t";
			
			my $MPOS = $1 if $line[120] =~ /MPOS\=(\d+)\;NALOD/;
			#my @MMQ = split /\,/, $MMQ;
			print OUT "$MPOS\t";
			#;NALOD=2.59;NLOD=113.77;POPAF=6.00;TLOD=128.92
			
			my $NALOD = $1 if $line[120] =~ /NALOD\=(\S+)\;NLOD/;
			#my @MMQ = split /\,/, $MMQ;
			if ($NALOD =~ /e/){
				print OUT "0\t";
			}
			else{
				print OUT "$NALOD\t";
			}
			
			
			my $NLOD = $1 if $line[120] =~ /NLOD\=(\S+)\;POPAF/;
			#my @MMQ = split /\,/, $MMQ;
			print OUT "$NLOD\t";
			
			my $POPAF = $1 if $line[120] =~ /POPAF\=(\S+)\;TLOD/;
			#my @MMQ = split /\,/, $MMQ;
			print OUT "$POPAF\t";
			
			my $TLOD = $1 if $line[120] =~ /TLOD\=(\S+)/;
			#my @MMQ = split /\,/, $MMQ;
			print OUT "$TLOD\t";
			
			if($line[5] eq "exonic"){
				print OUT "1\t0\t0\t0\t0\t0\t0\t";
			}
			elsif($line[5] eq "intronic"){
				print OUT "0\t1\t0\t0\t0\t0\t0\t";
			}
			elsif($line[5] =~ /UTR/){
				print OUT "0\t0\t1\t0\t0\t0\t0\t";
			}
			elsif($line[5] =~ /stream/){
				print OUT "0\t0\t0\t1\t0\t0\t0\t";
			}
			elsif($line[5] eq "intergenic"){
				print OUT "0\t0\t0\t0\t1\t0\t0\t";
			}
			elsif($line[5] eq "splicing"){
				print OUT "0\t0\t0\t0\t0\t1\t0\t";
			}
			else{
				print OUT "0\t0\t0\t0\t0\t0\t1\t";
			}
			my $size = keys %pat_total;
			if(exists $sig{$line[6]}{$line[1]} and exists $sig_total{$line[6]}){
				my $ratio1 = $sig{$line[6]}{$line[1]}/$size;
				my $ratio2 = $sig{$line[6]}{$line[1]}/$sig_total{$line[6]};
				print OUT "$sig{$line[6]}{$line[1]}\t";
				print OUT "$ratio1\t";
				print OUT "$ratio2\t";
			}
			else{
				print OUT "0\t";
				print OUT "1\t";
				print OUT "1\t";
			}
			
			if( (exists $wxs{$chr}{$pos} || exists $wxs{$chr}{$pos+1} || exists $wxs{$chr}{$pos+2} ||  exists $wxs{$chr}{$pos+3} || exists $wxs{$chr}{$pos+4} || exists $wxs{$chr}{$pos+5} || exists $wxs{$chr}{$pos-1} || exists $wxs{$chr}{$pos-2} || exists $wxs{$chr}{$pos-3} || exists $wxs{$chr}{$pos-4} || exists $wxs{$chr}{$pos-5} ) 
			&& (exists $wgs{$chr}{$pos} || exists $wgs{$chr}{$pos-1} || exists $wgs{$chr}{$pos-2} || exists $wgs{$chr}{$pos-3} || exists $wgs{$chr}{$pos-4} || exists $wgs{$chr}{$pos-5} || exists $wgs{$chr}{$pos+1} || exists $wgs{$chr}{$pos+2} || exists $wgs{$chr}{$pos+3} || exists $wgs{$chr}{$pos+4} || exists $wgs{$chr}{$pos+5} ) ){
				$share_count++;
				#print "$tmp{$chr}{$pos}\t";
				#print "wxs\twgs\t$rna{$chr}{$pos}\n";
				print OUT "1";
			}
			elsif(exists $wxs{$chr}{$pos} || exists $wxs{$chr}{$pos+1} || exists $wxs{$chr}{$pos+2} ||  exists $wxs{$chr}{$pos+3} || exists $wxs{$chr}{$pos+4} || exists $wxs{$chr}{$pos+5} || exists $wxs{$chr}{$pos-1} || exists $wxs{$chr}{$pos-2} || exists $wxs{$chr}{$pos-3} || exists $wxs{$chr}{$pos-4} || exists $wxs{$chr}{$pos-5}){
				$wxs_count++;
				#print "$tmp{$chr}{$pos}\t";
				#print "wxs\tno\t$rna{$chr}{$pos}\n";
				print OUT "1";
			}
			elsif(exists $wgs{$chr}{$pos} || exists $wgs{$chr}{$pos-1} || exists $wgs{$chr}{$pos-2} || exists $wgs{$chr}{$pos-3} || exists $wgs{$chr}{$pos-4} || exists $wgs{$chr}{$pos-5} || exists $wgs{$chr}{$pos+1} || exists $wgs{$chr}{$pos+2} || exists $wgs{$chr}{$pos+3} || exists $wgs{$chr}{$pos+4} || exists $wgs{$chr}{$pos+5}){
				$wgs_count++;
				#print "$tmp{$chr}{$pos}\t";
				#print "no\twgs\t$rna{$chr}{$pos}\n";
				print OUT "1";
			}
			else{
				#print "$tmp{$chr}{$pos}\t";
				#print "no\tno\t$rna{$chr}{$pos}\n";
				print OUT "0";
			}
			print OUT "\n";
		}
	}
	#print "$total_count\t$share_count\t$wxs_count\t$wgs_count\n";
	#print "$total_count\t$share_count\t$wxs_count\t$wgs_count\n";
}

