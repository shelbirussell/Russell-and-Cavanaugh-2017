use strict ;
use warnings ;
use Sort::Naturally ;
use Statistics::Histogram ;
use File::Basename ; 
use Data::Dumper ;
use Math::Counting ;

# Usage: perl intra-host_variants_from_pileup.pl Sv_sym_scaffold ref.fasta KB16gillMitoSym_stampy_realigned.pileup

my$scaff_keyword = $ARGV[0] ;
my$fasta = $ARGV[1] ;
my$pileup = $ARGV[2] ;

# error rate per bp
my$illumina_error_rate = 0.01 ;

print "scaffold keyword: ", $scaff_keyword, "\n" ;
print "reference fasta: ", $fasta, "\n" ;
print "pileup file: ", $pileup, "\n\n" ;

my$seqs = read_fasta($fasta, $scaff_keyword) ;
my%seqs = %{$seqs} ;

foreach my$scaff (keys %seqs) {
	print $scaff, "\n" ;
}

my($depths, $alleles, $min_alts) = parse_pileup($pileup, $scaff_keyword) ;
my@depths = @{$depths} ;
my%alleles = %{$alleles} ;
my@min_alts = @{$min_alts} ;

##########################################################################################
##############Filter by coverage depth and proximity to end of scaffold###################
### Calculate depth statistics to filter on
my$sum ;
$sum += $_ for @depths ;
my$mean = $sum / @depths ;

my$dev ;
foreach (@depths) {$dev += ($_ - $mean)**2 ;}
my$stdev = sqrt($dev / @depths) ;

my$stderr = $stdev / sqrt(@depths) ;

# 95% confidence interval bounds (Z=1.96); 99% confidence interval bounds (Z=2.576)
# Z * stdev or Z * stderr (stderr eliminates too many sites)
my$lowerbound = $mean - 1 * $stdev ;
my$upperbound = $mean + 1 * $stdev ;

# Minimum depth set to 100 manually
if ($mean-$stdev > 100 && $lowerbound < 100) {$lowerbound = 100 ;}

my$too_few_reads = 0 ;

my$out_of_bounds_count = 0 ;
my$end_of_scaff_count = 0 ;
my$counter = 0 ;

my@new_mins = () ;

foreach my$scaff (nsort keys %alleles) {
	foreach my$pos (nsort keys %{$alleles{$scaff}}) {
		if (!exists $alleles{$scaff}{$pos}{"DEPTH"}) {
			print $scaff, "\t", $pos, "\tdepth:", $alleles{$scaff}{$pos}{"DEPTH"}, "\n" ;
			delete $alleles{$scaff}{$pos} ;
			next ;
		}
		
		# exclude sites with coverage out of bounds
		if ($alleles{$scaff}{$pos}{"DEPTH"} > $upperbound || $alleles{$scaff}{$pos}{"DEPTH"} < $lowerbound) {
			$out_of_bounds_count ++ ;
			delete $alleles{$scaff}{$pos} ;
			next ;
		}
		
		#exclude sites within 100 bp of the end of a scaffold due to potential mapping error undetectable by coverage anomalies
		my$scaff_length = length($seqs{$scaff}) ;
		if ($pos > $scaff_length-100) {
			$end_of_scaff_count ++ ;
			delete $alleles{$scaff}{$pos} ;
			next ;
		}
		
		if ($pos < 100) {
			$end_of_scaff_count ++ ;
			delete $alleles{$scaff}{$pos} ;
			next ;
		}
			
		# skip site if 0 calls 
    	if ($alleles{$scaff}{$pos}{"REF_COUNT"} == 0 && $alleles{$scaff}{$pos}{"ALT_A_COUNT"} == 0 && $alleles{$scaff}{$pos}{"ALT_T_COUNT"} == 0 && $alleles{$scaff}{$pos}{"ALT_C_COUNT"} == 0 && $alleles{$scaff}{$pos}{"ALT_G_COUNT"} == 0 && $alleles{$scaff}{$pos}{"INDEL_COUNT"} == 0) {
			$too_few_reads ++ ;
			delete $alleles{$scaff}{$pos} ;
			next ;
		}	
	
		push @new_mins, $min_alts[$counter] ;
		$counter ++ ;
	
	}
}

my$min_sum = 0 ;
foreach my$min (@new_mins) {$min_sum += $min ;}
my$min_mean = $min_sum/@new_mins ;
	
my$min_dev = 0 ;
foreach (@new_mins) {$min_dev += ($_ - $min_mean)**2 ;}
my$min_stdev = sqrt($min_dev / @new_mins) ;

##########################################################################################
#########################Print biallelic variant site counts##############################
##########
my$AC_output = basename($pileup) ;
$AC_output =~ s/.pileup/${scaff_keyword}_all_AF_counts.txt/ ;
$AC_output =~ s/.filtered// ;
$AC_output =~ s/_stampy_realigned// ;

my$total_S = 0 ;
my$snps = 0 ;
my$indels = 0 ;

open OUT, ">$AC_output" or die "cannot open $AC_output\n" ;

foreach my$scaffold (nsort keys %alleles) {
	foreach my$position (sort {$a<=>$b} keys %{$alleles{$scaffold}}) {		
				
		my$ref_base = $alleles{$scaffold}{$position}{"REF_BASE"} ;
		my$ref_count = $alleles{$scaffold}{$position}{"REF_COUNT"} ;
		my%alts ;
		$alts{"A"} = $alleles{$scaffold}{$position}{"ALT_A_COUNT"} ;
		$alts{"T"} = $alleles{$scaffold}{$position}{"ALT_T_COUNT"} ;
		$alts{"C"} = $alleles{$scaffold}{$position}{"ALT_C_COUNT"} ;
		$alts{"G"} = $alleles{$scaffold}{$position}{"ALT_G_COUNT"} ;
		
		my$indel ;
		if ($alleles{$scaffold}{$position}{"INDEL"}) {
			$indel = $alleles{$scaffold}{$position}{"INDEL"} ;						
			$alts{$indel} = $alleles{$scaffold}{$position}{"INDEL_COUNT"} ;
			
			$total_S ++ ;
			$indels ++ ;
		}
		
		if ($ref_count > 0) {
			my@alt_bp = () ;
			foreach my$bp (keys %alts) {
				if ($alts{$bp} > 0) {push(@alt_bp, $bp) ;}				
			}
						
			if (scalar(@alt_bp) == 0) { next ;} # invariant site
			if (scalar(@alt_bp) > 1) { next ; print $scaffold, "\t", $position, "\t", "more than one alternate allele\n" ;}
			
			$total_S ++ ;
			$snps ++ ;
			
			print OUT $scaffold, "\t", $position, "\t", $ref_base, "\t", $ref_count, "\t", $alt_bp[0], "\t", $alts{$alt_bp[0]}, "\t", $alts{$alt_bp[0]}/($alts{$alt_bp[0]}+$ref_count), "\n" ;			
		}
		
		else {
			my@alt_bp ; 
			foreach my$bp (keys %alts) {
				if ($alts{$bp} > 0) {push(@alt_bp, $bp) ;}				
			}
			if (scalar(@alt_bp) < 2) { next ;} # invariant site
			if (scalar(@alt_bp) > 2) { next ; print $scaffold, "\t", $position, "\t", "more than one alternate allele\n" ;}

			$total_S ++ ;
			$snps ++ ;

			print OUT $scaffold, "\t", $position, "\t", $alt_bp[0], "\t", $alts{$alt_bp[0]}, "\t", $alt_bp[1], "\t", $alts{$alt_bp[1]}, "\t", $alts{$alt_bp[1]}/($alts{$alt_bp[1]}+$alts{$alt_bp[0]}), "\n" ;			
		}
	}
}

close OUT ;

print "Total intra-host segregating sites: ", $total_S, "\n" ;
print "Snps: ", $snps, "\n" ;
print "Indels: ", $indels, "\n" ;

##########################################################################################

##########################################################################################
#########################Print site by site multiallelic pi##############################
##### Iterate through reference sequence printing site-by-site pairwise diversity
my$pi_output = basename($pileup) ;
$pi_output =~ s/.pileup/${scaff_keyword}_pi.txt/ ;
$pi_output =~ s/.filtered// ;
$pi_output =~ s/_stampy_realigned// ;
open OUT, ">$pi_output" or die "cannot open $pi_output\n" ;

my$sum_pi = 0 ;

foreach my$scaff (nsort keys %seqs) {
	if ($scaff =~ m/$scaff_keyword/) {
		my$length = length($seqs{$scaff}) ;

		foreach my$pos (1..$length) {

			if ($alleles{$scaff}{$pos}) {
				my$ref = $alleles{$scaff}{$pos}{"REF_COUNT"} ;
				my$altA = $alleles{$scaff}{$pos}{"ALT_A_COUNT"} ;
				my$altT = $alleles{$scaff}{$pos}{"ALT_T_COUNT"} ;
				my$altC = $alleles{$scaff}{$pos}{"ALT_C_COUNT"} ;
				my$altG = $alleles{$scaff}{$pos}{"ALT_G_COUNT"} ;
				my$indel = $alleles{$scaff}{$pos}{"INDEL_COUNT"} ;
			
				# Multiallelic pi: over all site sum( over all alleles sum_i( Ji (n-Ji) ) / (n(n-1)) ) / L
				my@ACs = ($ref, $altA, $altT, $altC, $altG, $indel) ;
				my$AN = $ref + $altA + $altT + $altC + $altG + $indel ;
				my$sum_i = 0 ;
				foreach my$allele (@ACs) {$sum_i += ($allele * ($AN-$allele)) ;}
				my$site_pi = $sum_i / ($AN*($AN-1)) ;
			
				print OUT $scaff, "\t", $pos , "\t", $site_pi, "\n" ;
			
				$sum_pi += $site_pi ;
			}
		
			else {
				print OUT $scaff, "\t", $pos , "\t", "NA", "\n" ;
			}
		}
	}
}

close OUT ;

my$total_length = 0 ;
foreach my$scaffold (keys %seqs) {
	$total_length += length($seqs{$scaffold}) ;
}

print "Genome average pi: ", $sum_pi/$total_length, "\n\n" ;


##########################################################################################

print "Average minimum allele count required at a position: $min_mean +/- $min_stdev\n";
print "Sites with too few reads: $too_few_reads\n" ;
print "Coverage mean:", $mean, "\tstdev:", $stdev, "\tstderr:", $stderr, "\n" ;
print "Coverage cutoff lower_bound:", $lowerbound, "\tupper_bound:", $upperbound, "\n" ;
print "sites out of coverage bound: ", $out_of_bounds_count, "\n" ;
print "sites removed within 100 bp of scaffold ends: ", $end_of_scaff_count, "\n" ;
print "\n", get_histogram(\@depths) ;

sub parse_pileup {
	my$file = $_[0] ;
	my$keyword = $_[1] ;
	
	print "Reading pileup: ", $file, "\n\n" ;

	open PILEUP, "<$file" or die "cannot open $file\n\n" ;

	my@depths = () ;
	my%AC ;
	my@min_alts ;

	while (<PILEUP>) {
	
		chomp ; 
	
		if ($_ =~ m/^#/) {next;}
		
		my@split = split(/\t/, $_) ;
		
		if ($split[0] !~ m/$keyword/) {next ;}

		my$ref_base = $split[2] ;
		my$depth = $split[3] ;	

		### Set min depth to 1 (opposed to 0) to avoid issues with binomial calculation
		if ($depth <= 1) {next ;}

		push(@depths, $depth) ;

		#####
		## Call column must be processed in this order otherwise regex matches will be incorrect (e.g. all indel bp will be detected as snp calls)
		my$data = $split[4] ;
	
		# Remove mapping qualities (^ marks beginning of read and the ASCII of the character following it - 33 is is the mapping quality 
		$data =~ s/(\^.)//g ;		
		# Remove end of read symbol
		$data =~ s/\$//g ;
		
		# Count number of indels and remove from $data
		my$indel = 0 ;
		my@indels = () ;
		
		my%hash = () ;
		while ($data =~ /(\d+)/g) {
			$hash{$1} = 1 ;
			$indel ++ ;
		}		

		foreach my$k (keys %hash) {
			$data =~ s/[ATCGNatcgn\.,\*](-$k[ATCGNatcgn]{$k})//g ;		
			$data =~ s/[ATCGNatcgn\.,\*](\+$k[ATCGNatcgn]{$k})//g ;
			push(@indels, $1) ;
		}
						
		# Count reference bases deleted in deletions reported in previous lines (*)
		my$deleted = 0 ;
		while ($data =~ m/\*/gi) { $deleted ++ ;}
		
		# count number of reads with reference allele at position
		my$ref = 0 ;
		while ($data =~ m/[\.\,]/gi) { $ref ++ ; }
	
		# count number of other possible alternate alleles
		my$altA = 0;
		while ($data =~ m/[aA]/gi) {$altA ++ ;}

		my$altT = 0;
		while ($data =~ m/[tT]/gi) {$altT ++ ;}

		my$altC = 0;
		while ($data =~ m/[cC]/gi) {$altC ++ ;}

		my$altG = 0;
		while ($data =~ m/[gG]/gi) {$altG ++ ;}
		
		# If indel calls exist:
		if ($indel > 0) {
			
			# Capitalize indels so they match between strands
			$_ = uc foreach @indels ;
			
			# Check that all indels are the same, and if not, remove calls
			my%string = map{$_, 1} @indels ;
			if (keys %string == 1) {# all indels equal
			}

			else {$indel = 0 ;}
		}
	
		## Minimum minor allele count required based on errors at a position being binomially distributed
		## P(error) = (coverage choose error) * seq_error^error * (1-seq_error)^(coverage-error)
		my$min_alt ;
		my$prob_accurate = 0 ;
				
		foreach my$cov (0..$depth) {
			if ($prob_accurate >= 0.99) {last ;}
			
			$min_alt = $cov ;
			my$accurate = binom($depth, $cov) * ($illumina_error_rate)**$cov * (1-$illumina_error_rate)**($depth-$cov) ;
			$prob_accurate +=  $accurate ;

#			print $depth, "\t", $min_alt, "\t", $prob_accurate, "\n" ;
		}
	
		## Set absolute minimum to 5 
		if ($min_alt < 5) {$min_alt = 5 ;}
		
		push @min_alts, $min_alt ;
		
		## Output variant call if call > min_alt at site
		$AC{$split[0]}{$split[1]}{"REF_BASE"} = $ref_base ;
		$AC{$split[0]}{$split[1]}{"DEPTH"} = $depth ;

		if ($ref > $min_alt) {$AC{$split[0]}{$split[1]}{"REF_COUNT"} = $ref ;}
		else {$AC{$split[0]}{$split[1]}{"REF_COUNT"} = 0 ;}
		
		if ($altA > $min_alt) {$AC{$split[0]}{$split[1]}{"ALT_A_COUNT"} = $altA ;}
		else {$AC{$split[0]}{$split[1]}{"ALT_A_COUNT"} = 0 ;}
		
		if ($altT> $min_alt) {$AC{$split[0]}{$split[1]}{"ALT_T_COUNT"} = $altT ;}
		else {$AC{$split[0]}{$split[1]}{"ALT_T_COUNT"} = 0 ;}
		
		if ($altC > $min_alt){$AC{$split[0]}{$split[1]}{"ALT_C_COUNT"} = $altC ;}
		else {$AC{$split[0]}{$split[1]}{"ALT_C_COUNT"} = 0 ;}
		
		if ($altG > $min_alt) {$AC{$split[0]}{$split[1]}{"ALT_G_COUNT"} = $altG ;}
		else {$AC{$split[0]}{$split[1]}{"ALT_G_COUNT"} = 0 ;}

		if ($indel > $min_alt) {
			$AC{$split[0]}{$split[1]}{"INDEL_COUNT"} = $indel ;
			$AC{$split[0]}{$split[1]}{"INDEL"} = $indels[0] ;
		}
		
		else {$AC{$split[0]}{$split[1]}{"INDEL_COUNT"} = 0 ;}
		
	}
	
	close PILEUP ;

#	print "\nDone reading pileup\n\nSites with too few reads: $too_few_reads\nInvariant sites: $no_variant_count\n\n" ;

	## Delete snp calls overlapping or within 5 bp of an indel
	foreach my$scaffold (keys %AC) {
		my$positions = 0 ;
		my$deleted = 0 ;

		foreach my$position (keys %{$AC{$scaffold}}) {
			$positions ++ ;

			# check if position is an indel, and keep it if no snps were called there too
			if (exists $AC{$scaffold}{$position}{"INDEL"}) {
				my$alleles = 0 ;
				if ($AC{$scaffold}{$position}{"REF_COUNT"} > 0) {$alleles ++ ;}
				if ($AC{$scaffold}{$position}{"ALT_A_COUNT"} > 0) {$alleles ++ ;}
				if ($AC{$scaffold}{$position}{"ALT_T_COUNT"} > 0) {$alleles ++ ;}
				if ($AC{$scaffold}{$position}{"ALT_C_COUNT"} > 0) {$alleles ++ ;}
				if ($AC{$scaffold}{$position}{"ALT_G_COUNT"} > 0) {$alleles ++ ;}
				
				if ($alleles >1) {delete $AC{$scaffold}{$position} ;}
				else {next ;}
			}

			else {
				my@range = ($position-5..$position+5) ;
				my$to_delete ;
				foreach my$i (@range) {
					# delete snp site if within 5 bp of indel (low and high confidence indels)
					if ($AC{$scaffold}{$position}{"INDEL"}) {					
						$to_delete = "yes" ;
					}
				}
				if ($to_delete && $to_delete eq "yes") {
					$deleted ++ ;
					delete $AC{$scaffold}{$position} ;
				}
			}
		}
		print "Positions deleted on $scaffold within 5bp of indel: $deleted\n" ;
	}
	
	return (\@depths, \%AC, \@min_alts) ;
}

sub read_fasta {
	print "Reading fasta: ", $_[0], "\n\n" ;

    open FASTA, "<$_[0]" or die "can not read fasta $_[0]\n" ;

    my$keyword = $_[1] ;
    my%seqs ;
    my$header ;
    my$seq ;
	my$count = 0 ;

    while (<FASTA>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
                $seqs{$header} = $seq ;
            }

            chomp ;

            my$header_line = $_ ;

			my@header_elements = split(/\t/, $header_line) ;
            $header = $header_elements[0] ;
            $header =~ s/^>// ;
            $header =~ s/\s+$// ;

			$count ++ ;
			
            $seq = "" ;
        }

        else {
            $_ =~ s/\s+//g ;
            $seq .= $_ ;
        }
    }

    close FASTA ;

    if ($seq) {
        $seqs{$header} = $seq ;
    }

    #delete non-keyword scaffolds
    foreach my$scaff (keys %seqs) {
    	if ($scaff !~ m/$keyword/) {
    		delete $seqs{$scaff} ;
    	}
    }
    
    return \%seqs ;
}

sub binom {
    my( $n, $r ) = @_;
    return unless defined $n && $n =~ /^\d+$/ && defined $r && $r =~ /^\d+$/;
    my $product = 1;
    while( $r > 0 ) {
        $product *= $n--;
        $product /= $r--;
    }
    return $product;
}
