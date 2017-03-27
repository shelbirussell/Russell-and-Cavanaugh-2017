use strict ; 
use warnings ;
use File::Basename ;

my$window = $ARGV[0] ;
my$file = $ARGV[1] ;
my$out = basename($file) ;
$out =~ s/.txt/${window}.txt/ ;

open FILE, "<$file" or die "cannot open $file\n" ;
open OUT, ">$out" or die "cannot open $out\n" ;

my%lengths ;

#$lengths{"Sv_sym_scaffold1"} = 1213831 ;
#$lengths{"Sv_sym_scaffold2"} = 892555 ;
#$lengths{"Sv_sym_scaffold3"} = 537613 ;
#$lengths{"Sv_sym_scaffold4"} = 28016 ;

my$window_counter = 0 ;
my$div_counter = 0 ;
my$sum = 0 ;
my$snps = 0 ;
my$prev_scaff ;
my$prev_pos ;

while (<FILE>) {
    if ($_ =~ m/^#/) { next ;}
    chomp ;
	my@split = split(/\t/, $_) ;
		
    if (! $prev_scaff) {$prev_scaff = $split[0] ; $prev_pos = $split[1] - 1 ;}

	if ($split[0] ne $prev_scaff) {
        print OUT $prev_scaff, "\t", $prev_pos, "\t", $div_counter/$window_counter, "\t", $snps, "\t", $sum/$window_counter, "\n" ;
        
		$window_counter = 1 ;
        $div_counter = 1 ;

        if ($split[2] =~ m/\d+/) {
            $sum = $split[2] ;
	        if ($split[2] > 0) {$snps = 1 ;}
        }
        
		else {
			$sum = 0 ;
			$snps = 0 ;
		}

		$prev_scaff = $split[0] ;
        $prev_pos = $split[1] ;
	}
	
	elsif ($window_counter == $window) {
		print OUT $split[0], "\t", $split[1], "\t", $div_counter/$window_counter, "\t", $snps, "\t", $sum/$window_counter, "\n" ;
		$window_counter = 1 ;
        $div_counter = 1 ;
        $prev_pos = $split[1] ;
        
        if ($split[2] =~ m/\d+/) {
            $sum = $split[2] ;
            if ($split[2] > 0) {$snps = 1 ;}
        }
        
		else {
			$sum = 0 ;
			$snps = 0 ;
		}
	}
	
    
    else {
        if ($split[2] =~ m/\d+/) {
            $div_counter ++ ;
            $sum += $split[2] ;
            if ($split[2] > 0) {$snps ++ ;}
        }

        $window_counter ++ ;
        $prev_pos = $split[1] ;
	}
}

close FILE ;

if ($window_counter > 0) {
    print OUT $prev_scaff, "\t", $prev_pos, "\t", $div_counter/$window_counter, "\t", $snps, "\t", $sum/$window_counter, "\n" ;
}
    
close OUT ;
