use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

#Usage: perl parse_intrahost_sites.pl Sv_biallelic_AF.txt KB16_candidate_sites.txt
# Use this to plot sites in particular samples (e.g. "sites.txt") along the genome (e.g. "subpop_biallelic_AF.txt")
# Also useful for marking haplotypes supported by reads

my$subpop_sites = $ARGV[0] ;
my$selected_sites = $ARGV[1] ;

my$output = basename($selected_sites) ;
$output =~ s/MitoSym.filtered(_\d+\.\d+_percent_)AF.txt/$1plot_AF.txt/ ;
$output =~ s/read_haplotypes.txt/haplotypes_plot_AF.txt/ ;

open SUBPOP, "<$subpop_sites" or die "cannot open $subpop_sites\n" ;

my%sites ; 

while (<SUBPOP>) {
	chomp $_ ;
	
	next if ($_ =~ m/^#/) ;
	
	my@split = split(/\t/, $_) ;
	
	if ($split[4]) {
		if ($split[4] == 0.00 || $split[4] == 1.00) {next ;}
	}
	
	$sites{$split[0]}{$split[1]} = "NA" ;
	
}

close SUBPOP ;

open INTRAHOST, "<$selected_sites" or die "cannot open $selected_sites\n" ;

while (<INTRAHOST>) {
	chomp ;
	
	next if ($_ =~ m/^#/) ;
	
	my@split = split(/\t/, $_) ;
	
	$sites{$split[0]}{$split[1]} = 1 ;
		
}	

close INTRAHOST ;

open OUT, ">$output" ;

foreach my$scaff (nsort keys %sites) {
	foreach my $pos (nsort keys %{$sites{$scaff}}) {
		print OUT $scaff, "\t", $pos, "\t", $sites{$scaff}{$pos}, "\n" ;
	}
}

close OUT ;
