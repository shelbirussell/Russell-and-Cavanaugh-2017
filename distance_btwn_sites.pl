use strict ;
use warnings ;

my$positions = $ARGV[0] ;

open POSITIONS, "<$positions" ;

my$above_350bp = 1 ;
my$below_350bp = 1 ;
my$total = 0 ;

my$scaffold ;
my$pos ;

while (<POSITIONS>) {
	chomp ;
	my@split = split(/\t/, $_) ;	
	
	$total ++ ;

	if (!$scaffold) {
		$scaffold = $split[0] ;
		$pos = $split[1] ;
		next ;
	}
	
	if ($split[0] eq $scaffold) {
		my$distance = $split[1] - $pos ;
		
		if ($distance > 350) {$above_350bp ++ ;}
		else {$below_350bp ++ ;}
	
		$scaffold = $split[0] ;
		$pos = $split[1] ;	
	}

	else {
		$scaffold = $split[0] ;
		$pos = $split[1] ;
	}	
}

close POSITIONS ;

print "Pairwise positions greater than 350 bp apart: ", $above_350bp, "\n" ;
print "Pairwise positions less than 350 bp apart: ", $below_350bp, "\n" ;
print "Total positions: ", $total, "\n" ;
print "Fraction of pairwise positions less than 350 bp apart: ", $below_350bp/$total, "\n" ;