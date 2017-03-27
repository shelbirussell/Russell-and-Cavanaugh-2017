use strict ;
use warnings ;
use Sort::Naturally ;
#use Statistics::PointEstimation ;
use Statistics::Histogram ;
use File::Basename; 

# Usage: perl allele_counts-segsites.pl SvKB_biallelic_AF.txt KB16gillMitoSym_all_AF_counts.txt

my$AF_file = $ARGV[0] ;
my$AF_counts = $ARGV[1] ;
my$output_pop = basename($AF_counts) ;
$output_pop =~ s/.txt/_segsites.txt/ ;
my$output_subpop = basename($AF_counts) ;
$output_subpop =~ s/.txt/_subpopsegsites.txt/ ;

print basename($AF_counts), "\n" ;

my$invariant_in_pop = 0 ;
my$absent_in_subpop = 0 ;
my$different_alleles = 0 ;

open AF, "<$AF_file" or die "cannot open $AF_file\n" ;

my%AF ;

while (<AF>) {
	chomp ;

	my@split = split(/\t/, $_) ;
	
	$AF{$split[0]}{$split[1]}{"REF"} = $split[2] ;
	$AF{$split[0]}{$split[1]}{"ALT"} = $split[3] ;
	$AF{$split[0]}{$split[1]}{"AF"} = $split[4] ;
}

close AF ;

open ALLELES, "<$AF_counts" or die "cannot open $AF_counts\n" ;

while (<ALLELES>) {
	chomp ; 
	
	my@split = split(/\t/, $_) ;
	my$scaffold = $split[0] ;
	my$position = $split[1] ;
	my$base1 = $split[2] ;
	my$count1 = $split[3] ;
	my$base2 = $split[4] ;
	my$count2 = $split[5] ;
	
	if (exists $AF{$scaffold}{$position}) {
#		print substr($AF{$scaffold}{$position}{"REF"},1), "\t", $base1, "\t", substr($AF{$scaffold}{$position}{"ALT"}, 1), "\t", $base2, "\n" ;

		# If + --> insertion in alternate
		if ($base2 =~ m/\+[0-9]+([ATCGatcg][ATCGatcg]+)/) {
			if ($1 eq substr($AF{$scaffold}{$position}{"ALT"}, 1) && $base1 eq $AF{$scaffold}{$position}{"REF"}) {
				$AF{$split[0]}{$split[1]}{"REF_COUNT"} = $count1 ;
				$AF{$split[0]}{$split[1]}{"ALT_COUNT"} = $count2 ;		
			}
		}	

		# If - --> deletion in reference
		elsif ($base1 =~ m/-[0-9]+([ATCGatcg][ATCGatcg]+)/) {
			if ($1 eq substr($AF{$scaffold}{$position}{"REF"}, 1) && $base2 eq $AF{$scaffold}{$position}{"ALT"}) {
				$AF{$split[0]}{$split[1]}{"REF_COUNT"} = $count1 ;
				$AF{$split[0]}{$split[1]}{"ALT_COUNT"} = $count2 ;		
			}
		}
		elsif ($base2 =~ m/-[0-9]+([ATCGatcg][ATCGatcg]+)/) {
			if ($1 eq substr($AF{$scaffold}{$position}{"REF"}, 1) && $base1 eq $AF{$scaffold}{$position}{"ALT"}) {
#				print $1, "\t", $AF{$scaffold}{$position}{"REF"}, "\t", substr($AF{$scaffold}{$position}{"REF"}, 1), "\n" ;
				$AF{$split[0]}{$split[1]}{"REF_COUNT"} = $count1 ;
				$AF{$split[0]}{$split[1]}{"ALT_COUNT"} = $count2 ;		
			}
		}
					
		elsif ($base1 eq $AF{$scaffold}{$position}{"REF"} && $base2 eq $AF{$scaffold}{$position}{"ALT"}) {
			$AF{$split[0]}{$split[1]}{"REF_COUNT"} = $count1 ;
			$AF{$split[0]}{$split[1]}{"ALT_COUNT"} = $count2 ;
		}
		
		elsif ($base1 eq $AF{$scaffold}{$position}{"ALT"} && $base2 eq $AF{$scaffold}{$position}{"REF"}) {
			$AF{$split[0]}{$split[1]}{"REF_COUNT"} = $count1 ;
			$AF{$split[0]}{$split[1]}{"ALT_COUNT"} = $count2 ;
		}

		else {
			$different_alleles ++ ;
		}
	}
}

foreach my$scaffold (keys %AF) {
	foreach my$position (keys %{$AF{$scaffold}}) {
#		print $scaffold, "\t", $position, "\t", $AF{$scaffold}{$position}{"REF_COUNT"}, "\t", $AF{$scaffold}{$position}{"ALT_COUNT"}, "\t", $AF{$scaffold}{$position}{"AF"}, "\n" ;
	}
}

close ALLELES ;

my$segS = 0 ;
my$subpopS = 0 ;
my$segsnps = 0 ;
my$subpopsnps = 0 ;
my$segindels = 0 ;
my$subpopindels = 0 ;

open OUT_POP, ">$output_pop" or die "cannot open $output_pop\n" ;
open OUT_SUBPOP, ">$output_subpop" or die "cannot open $output_subpop\n" ;

foreach my$scaff (nsort keys %AF) {
	foreach my$pos (sort {$a<=>$b} keys %{$AF{$scaff}}) {		
		
		## Skip site if not called in sample
		if (!exists $AF{$scaff}{$pos}{"REF_COUNT"} && !exists $AF{$scaff}{$pos}{"ALT_COUNT"}) {
			next ;
		}
		
		## Exclude sites that aren't variant in this intra-host population (0,0,0)
		if ($AF{$scaff}{$pos}{"REF_COUNT"} == 0 && $AF{$scaff}{$pos}{"ALT_COUNT"} == 0) {
			next ;
		}
		
		if ($AF{$scaff}{$pos}{"REF"} =~ m/[ATCGatcg]{2,}/ || $AF{$scaff}{$pos}{"ALT"} =~ m/[ATCGatcg]{2,}/) {
			$segS ++ ;
			$segindels ++ ;
		}

		else {
			$segS ++ ;
			$segsnps ++ ;
		}
		
		print OUT_POP $scaff, "\t", $pos , "\t", $AF{$scaff}{$pos}{"REF"}, "\t", $AF{$scaff}{$pos}{"REF_COUNT"}, "\t", $AF{$scaff}{$pos}{"ALT"}, "\t", $AF{$scaff}{$pos}{"ALT_COUNT"}, "\t", $AF{$scaff}{$pos}{"ALT_COUNT"}/($AF{$scaff}{$pos}{"ALT_COUNT"}+$AF{$scaff}{$pos}{"REF_COUNT"}), "\t", $AF{$scaff}{$pos}{"AF"}, "\n" ;


		## Exclude sites from subpop file that aren't in the subpopulation
		if ($AF{$scaff}{$pos}{"AF"} eq "0.00") {
			$absent_in_subpop ++ ;
			next ;
		}
		
		if ($AF{$scaff}{$pos}{"REF"} =~ m/[ATCGatcg]{2,}/ || $AF{$scaff}{$pos}{"ALT"} =~ m/[ATCGatcg]{2,}/) {
			$subpopS ++ ;
			$subpopindels ++ ;
		}

		else {
			$subpopS ++ ;
			$subpopsnps ++ ;
		}

		print OUT_SUBPOP $scaff, "\t", $pos , "\t", $AF{$scaff}{$pos}{"REF"}, "\t", $AF{$scaff}{$pos}{"REF_COUNT"}, "\t", $AF{$scaff}{$pos}{"ALT"}, "\t", $AF{$scaff}{$pos}{"ALT_COUNT"}, "\t", $AF{$scaff}{$pos}{"AF"}, "\n" ;
	}
}

close OUT_POP ;
close OUT_SUBPOP ;

print "sites with different alleles in intra-host population than between-host populations: ", $different_alleles, "\n" ;
print "sites fixed in subpopulation:", $invariant_in_pop, "\n" ;
print "sites absent in subpopulation:", $absent_in_subpop, "\n" ;

print "Total segregating sites in population: ", $segS, "\n" ;
print "snps in population: ", $segsnps, "\n" ;
print "indels in population: ", $segindels, "\n" ;
print "Total segregating sites in subpopulation: ", $subpopS, "\n" ;
print "snps in subpopulation: ", $subpopsnps, "\n" ;
print "indels in subpopulation: ", $subpopindels, "\n" ;
