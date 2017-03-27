use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# usage: perl get_alleles_by_AF.pl AF_counts.txt 0.20 0.50

my$AF_counts = $ARGV[0] ;
my$output_prefix = basename($AF_counts) ;
$output_prefix =~ s/all_// ;
$output_prefix =~ s/AF_counts_// ;
$output_prefix =~ s/.txt// ;

my@output_freq = () ;
foreach my$i (1..$#ARGV) {push @output_freq, $ARGV[$i] ;}
print join(",", @output_freq), "\n" ;

open COUNTS, "<$AF_counts" or die "cannot open $AF_counts\n" ;

my%folded_afs ;

while (<COUNTS>) {
	chomp ;
	if ($_ =~ m/^#/) {next;}
	my@split = split(/\t/, $_) ;
	
	my$scaffold = $split[0] ;
	my$pos = $split[1] ;
	my$baseA = $split[2] ;
	my$CA = $split[3] ;
	my$basea = $split[4] ;
	my$Ca = $split[5] ;
	
	my$af = $Ca / ($Ca + $CA) ;

	#fold site frequency
	if ($af > 0.5) {$af = 1 - $af ;}
	
	$folded_afs{$scaffold}{$pos}{"AF"} = $af ;	
	$folded_afs{$scaffold}{$pos}{"REF"} = $baseA ;	
	$folded_afs{$scaffold}{$pos}{"REF_count"} = $CA ;	
	$folded_afs{$scaffold}{$pos}{"ALT"} = $basea ;	
	$folded_afs{$scaffold}{$pos}{"ALT_count"} = $Ca ;	
}

foreach my$freq (@output_freq) {
	open OUT, ">${output_prefix}${freq}_AF.txt" ;
	
	foreach my$scaff (nsort keys %folded_afs) {
		foreach my$pos (nsort keys %{$folded_afs{$scaff}}) {
			my$max = $freq + 0.04 ;
			my$min = $freq - 0.04 ;
#			print "max:", $max, " min:", $min, "\n" ;
			if ($folded_afs{$scaff}{$pos}{"AF"} > $min && $folded_afs{$scaff}{$pos}{"AF"} < $max) {
				print OUT $scaff, "\t", $pos, "\t", $folded_afs{$scaff}{$pos}{"REF"}, "\t", $folded_afs{$scaff}{$pos}{"REF_count"}, "\t", $folded_afs{$scaff}{$pos}{"ALT"}, "\t", $folded_afs{$scaff}{$pos}{"ALT_count"}, "\t", $folded_afs{$scaff}{$pos}{"AF"}, "\n" ;
			}
		}
	}
	close OUT ;
}
