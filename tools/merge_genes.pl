use strict;
use warnings;

my $gtf = $ARGV[0];
my $out_files = $ARGV[1];

open(IN, "$gtf") or die "Cannot open $gtf file for reading: $!\n";

my %data;
while (<IN>) {
    next if $_ =~ /^#/;	
    $_ =~ s/\s+$//;
    my @line = split /\t/, $_;
	next unless $line[2] eq 'gene'; 

    push @{ $data{$line[0]} }, [$line[3], $line[4]];
}

close(IN);

open(OUT, ">$out_files") or die "Cannot open $out_files file for writing: $!\n";

for my $chr (sort keys %data) {
    my @ranges = sort { $a->[0] <=> $b->[0] } @{ $data{$chr} };
    my @merged = ();

    for my $range (@ranges) {
        if (!@merged or $merged[-1]->[1] < $range->[0]) {
            push @merged, $range;
        }
        else {
            $merged[-1]->[1] = $range->[1] if $range->[1] > $merged[-1]->[1];
        }
    }

    for my $range (@merged) {
        print OUT "$chr:$range->[0]-$range->[1]\n";
    }
}

close(OUT);