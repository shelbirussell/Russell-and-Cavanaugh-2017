use strict ;
use warnings ;
use Sort::Naturally ;
use File::Basename ;

# Usage: perl intra-host_variants_from_mpileup.pl Sv_sym_scaffold ref.fasta KB16gillMitoSym_stampy_realigned.pileup

my$scaff_keyword = $ARGV[0] ;
my$fasta = $ARGV[1] ;
my$pileup_file = $ARGV[2] ;

# error rate per bp
my$illumina_error_rate = 0.01 ;

print "scaffold keyword: ", $scaff_keyword, "\n" ;
print "reference fasta: ", $fasta, "\n" ;
print "pileup file: ", $pileup_file, "\n\n" ;

my$seqs = read_fasta($fasta, $scaff_keyword) ;
my%seqs = %{$seqs} ;
print "FASTA\n" ;
foreach my$scaffold (keys %seqs) {
	print $scaffold, "\t", length($seqs{$scaffold}), "\n" ;
}

my$alleles = parse_pileup($pileup_file, $scaff_keyword) ;
my%alleles = %{$alleles} ;


## Print statistics to an organized table for easier comparison/viewing
my$stat_table = basename($pileup_file) ;
$stat_table =~ s/.pileup/${scaff_keyword}_pop_stats.txt/ ;

open SUMMARY, ">$stat_table" or die "cannot open $stat_table\n" ;
print SUMMARY "sample\ttotalS\tSNPs\tindels\taveragePi\n" ;


################################################################################################################################
#############Filter each sample's sites around indels, by coverage depth, and by proximity to end of scaffold###################
foreach my$sample (nsort keys %alleles) {
	my$total_sites = 0 ;
	my@depths = () ;

	foreach my$scaffold (nsort keys %{$alleles{$sample}}) {
		my$deleted = 0 ;

		foreach my$position (nsort keys %{$alleles{$sample}{$scaffold}}) {
			## Delete snp calls overlapping
			# check if position is an indel, and keep it if no snps were called there too
			if (exists $alleles{$sample}{$scaffold}{$position}{"INDEL"}) {
				my$al = 0 ;
				if ($alleles{$sample}{$scaffold}{$position}{"REF_COUNT"} > 0) {$al ++ ;}
				if ($alleles{$sample}{$scaffold}{$position}{"ALT_A_COUNT"} > 0) {$al ++ ;}
				if ($alleles{$sample}{$scaffold}{$position}{"ALT_T_COUNT"} > 0) {$al ++ ;}
				if ($alleles{$sample}{$scaffold}{$position}{"ALT_C_COUNT"} > 0) {$al ++ ;}
				if ($alleles{$sample}{$scaffold}{$position}{"ALT_G_COUNT"} > 0) {$al ++ ;}

				if ($al >1) {delete $alleles{$sample}{$scaffold}{$position} ;}
				else {
					push @depths, $alleles{$sample}{$scaffold}{$position}{"DEPTH"} ;
					$total_sites ++ ;
				}
			}

			## Delete snps within 5 bp of an indel
			else {
				my@range = ($position-5..$position+5) ;
				my$to_delete ;
				foreach my$i (@range) {
				# delete snp site if within 5 bp of indel (low and high confidence indels)
					if (exists $alleles{$sample}{$scaffold}{$i} && $alleles{$sample}{$scaffold}{$i}{"INDEL"}) {$to_delete = "yes" ;}
				}
				if ($to_delete && $to_delete eq "yes") {
					$deleted ++ ;
					delete $alleles{$sample}{$scaffold}{$position} ;
				}

					else {
						push @depths, $alleles{$sample}{$scaffold}{$position}{"DEPTH"} ;
						$total_sites ++ ;
					}
				}
			}
			print "Positions deleted in sample #$sample on $scaffold within 5bp of indel: $deleted\n" ;
		}

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

	my@min_alts = () ;

	foreach my$scaff (nsort keys %{$alleles{$sample}}) {
		foreach my$pos (nsort keys %{$alleles{$sample}{$scaff}}) {

			# exclude sites with coverage out of bounds
			if (!exists $alleles{$sample}{$scaff}{$pos}{"DEPTH"} || $alleles{$sample}{$scaff}{$pos}{"DEPTH"} > $upperbound || $alleles{$sample}{$scaff}{$pos}{"DEPTH"} < $lowerbound) {
				$out_of_bounds_count ++ ;
				delete $alleles{$sample}{$scaff}{$pos} ;
				next ;
			}

			#exclude sites within 100 bp of the end of a scaffold due to potential mapping error undetectable by coverage anomalies
			my$scaff_length = length($seqs{$scaff}) ;
			if ($pos > $scaff_length-100) {
				$end_of_scaff_count ++ ;
				delete $alleles{$sample}{$scaff}{$pos} ;
				next ;
			}

			if ($pos < 100) {
				$end_of_scaff_count ++ ;
				delete $alleles{$sample}{$scaff}{$pos} ;
				next ;
			}

			# skip site if 0 calls
	    	if ($alleles{$sample}{$scaff}{$pos}{"REF_COUNT"} == 0 && $alleles{$sample}{$scaff}{$pos}{"ALT_A_COUNT"} == 0 && $alleles{$sample}{$scaff}{$pos}{"ALT_T_COUNT"} == 0 && $alleles{$sample}{$scaff}{$pos}{"ALT_C_COUNT"} == 0 && $alleles{$sample}{$scaff}{$pos}{"ALT_G_COUNT"} == 0 && $alleles{$sample}{$scaff}{$pos}{"INDEL_COUNT"} == 0) {
				$too_few_reads ++ ;
				delete $alleles{$sample}{$scaff}{$pos} ;
				next ;
			}

			push @min_alts, $alleles{$sample}{$scaff}{$pos}{"MIN_ALT"} ;
		}
	}

	my$min_sum = 0 ;
	foreach my$min (@min_alts) {$min_sum += $min ;}

#	##skip to next sample if @min_alts is empty
	if (@min_alts == 0) {next ;}

	my$min_mean = $min_sum/(scalar@min_alts) ;

	my$min_dev = 0 ;
	foreach (@min_alts) {$min_dev += ($_ - $min_mean)**2 ;}
	my$min_stdev = sqrt($min_dev / (scalar@min_alts)) ;

	##########################################################################################
	#########################Print biallelic variant site counts##############################
	##########
	my$AC_output = basename($pileup_file) ;
	$AC_output =~ s/.pileup/${sample}_${scaff_keyword}_all_AF_counts.txt/ ;

	open OUT, ">$AC_output" or die "cannot open $AC_output\n" ;

	my$total_S = 0 ;
	my$snps = 0 ;
	my$indels = 0 ;


	foreach my$scaffold (nsort keys %{$alleles{$sample}}) {
		foreach my$position (sort {$a<=>$b} keys %{$alleles{$sample}{$scaffold}}) {

			## Get allele and counts at position
			my$ref_base = $alleles{$sample}{$scaffold}{$position}{"REF_BASE"} ;
			my$ref_count = $alleles{$sample}{$scaffold}{$position}{"REF_COUNT"} ;
			my$indel ;
			my%alts ;
			$alts{"A"} = $alleles{$sample}{$scaffold}{$position}{"ALT_A_COUNT"} ;
			$alts{"T"} = $alleles{$sample}{$scaffold}{$position}{"ALT_T_COUNT"} ;
			$alts{"C"} = $alleles{$sample}{$scaffold}{$position}{"ALT_C_COUNT"} ;
			$alts{"G"} = $alleles{$sample}{$scaffold}{$position}{"ALT_G_COUNT"} ;

			if ($alleles{$sample}{$scaffold}{$position}{"INDEL"}) {
				$indel = $alleles{$sample}{$scaffold}{$position}{"INDEL"} ;
				$alts{$indel} = $alleles{$sample}{$scaffold}{$position}{"INDEL_COUNT"} ;
			}

			## Count segregating alleles at position
			my@alt_bp = () ;
			foreach my$bp (keys %alts) {if ($alts{$bp} > 0) {push(@alt_bp, $bp) ;}}

			## if one of the alleles is the reference
			if ($ref_count > 0) {
				if (scalar(@alt_bp) == 0) { next ;} # invariant site
				elsif (scalar(@alt_bp) > 1) {print $sample, "\t", $scaffold, "\t", $position, "\t", "more than one alternate allele\n" ; next ;}
				else {
					print OUT $scaffold, "\t", $position, "\t", $ref_base, "\t", $ref_count, "\t", $alt_bp[0], "\t", $alts{$alt_bp[0]}, "\t", $alts{$alt_bp[0]}/($alts{$alt_bp[0]}+$ref_count), "\n" ;

					## since site is biallelic and encodes the reference base, the other site can either be an indel or a snp
					if ($alleles{$sample}{$scaffold}{$position}{"INDEL"}) {$indels ++ ;}
					else {$snps ++ ;}
					$total_S ++ ;
				}
			}

			## if none of the alleles are the reference
			else {
				if (scalar(@alt_bp) < 2) { next ;} # invariant site
				elsif (scalar(@alt_bp) > 2) {print $sample, "\t", $scaffold, "\t", $position, "\t", "more than one alternate allele\n" ; next ;}
				else {
					print OUT $scaffold, "\t", $position, "\t", $alt_bp[0], "\t", $alts{$alt_bp[0]}, "\t", $alt_bp[1], "\t", $alts{$alt_bp[1]}, "\t", $alts{$alt_bp[1]}/($alts{$alt_bp[1]}+$alts{$alt_bp[0]}), "\n" ;

					## biallelic site without the reference base could encode two snps or a snp and an indel
					## we don't want to count the site twice, so the indel will take priority
					if ($alleles{$sample}{$scaffold}{$position}{"INDEL"}) {$indels ++ ;}
					else {$snps ++ ;}
					$total_S ++ ;
				}
			}
		}
	}

	close OUT ;

	print "\nSAMPLE: ", $sample, "\n" ;
	print "Total intra-host segregating sites: ", $total_S, "\n" ;
	print "Snps: ", $snps, "\n" ;
	print "Indels: ", $indels, "\n" ;

	##########################################################################################

	##########################################################################################
	#########################Print site by site multiallelic pi##############################
	##### Iterate through reference sequence printing site-by-site pairwise diversity
	my$pi_output = basename($pileup_file) ;
	$pi_output =~ s/.pileup/${sample}_${scaff_keyword}_pi.txt/ ;

	open OUT, ">$pi_output" or die "cannot open $pi_output\n" ;

	my$sum_pi = 0 ;

	foreach my$scaff (nsort keys %seqs) {
		if ($scaff =~ m/$scaff_keyword/) {
			my$length = length($seqs{$scaff}) ;

			foreach my$pos (1..$length) {

				if ($alleles{$sample}{$scaff}{$pos}) {
					my$ref = $alleles{$sample}{$scaff}{$pos}{"REF_COUNT"} ;
					my$altA = $alleles{$sample}{$scaff}{$pos}{"ALT_A_COUNT"} ;
					my$altT = $alleles{$sample}{$scaff}{$pos}{"ALT_T_COUNT"} ;
					my$altC = $alleles{$sample}{$scaff}{$pos}{"ALT_C_COUNT"} ;
					my$altG = $alleles{$sample}{$scaff}{$pos}{"ALT_G_COUNT"} ;
					my$indel = $alleles{$sample}{$scaff}{$pos}{"INDEL_COUNT"} ;

					# Multiallelic pi: over all site sum( over all alleles sum_i( Ji (n-Ji) ) / (n(n-1)) ) / L
					my@ACs = ($ref, $altA, $altT, $altC, $altG, $indel) ;
					my$AN = $ref + $altA + $altT + $altC + $altG + $indel ;
					my$sum_i = 0 ;
					foreach my$allele (@ACs) {$sum_i += ($allele * ($AN-$allele)) ;}
					my$site_pi = $sum_i / ($AN*($AN-1)) ;

					print OUT $scaff, "\t", $pos , "\t", $site_pi, "\n" ;

					$sum_pi += $site_pi ;
				}

				else {print OUT $scaff, "\t", $pos , "\t", "NA", "\n" ;}
			}
		}
	}

	close OUT ;

	my$total_length = 0 ;
	foreach my$scaffold (keys %seqs) {$total_length += length($seqs{$scaffold}) ;}

	print "SAMPLE: ", $sample, "\n" ;
	print "Genome average pi: ", $sum_pi/$total_length, "\n" ;


	##########################################################################################
	print "SAMPLE: ", $sample, "\n" ;
	print "Average minimum allele count required at a position: $min_mean +/- $min_stdev\n";
	print "Sites with too few reads: $too_few_reads\n" ;
	print "Coverage mean:", $mean, "\tstdev:", $stdev, "\tstderr:", $stderr, "\n" ;
	print "Coverage cutoff lower_bound:", $lowerbound, "\tupper_bound:", $upperbound, "\n" ;
	print "Total number of sites: ", $total_sites, "\n" ;
	print "sites out of coverage bound: ", $out_of_bounds_count, "\n" ;
	print "sites removed within 100 bp of scaffold ends: ", $end_of_scaff_count, "\n" ;
#	print "\n", get_histogram(\@depths) ;
	print "\n\n" ;

	########################################################################################
	#### Print table of S and pi values
	print SUMMARY $sample, "\t", $total_S, "\t", $snps, "\t", $indels, "\t", $sum_pi/$total_length, "\n" ;

}

close SUMMARY ;



sub parse_pileup {
	my$file = $_[0] ;
	my$keyword = $_[1] ;

	print "Reading pileup: ", $file, "\n\n" ;

	open PILEUP, "<$file" or die "cannot open $file\n\n" ;

	my%mpileup ;

	while (<PILEUP>) {

		chomp ;

		if ($_ =~ m/^#/) {next;}

		my@split = split(/\t/, $_) ;

		if ($split[0] !~ m/$keyword/) {next ;}

		my$scaffold = $split[0] ;
		my$position = $split[1] ;
		my$ref_base = $split[2] ;

		##mpileups (multi sample pileups) have 3 columns for each sample, starting at column #3.
		##for each sample, these 3 columns contain: coverage depth, mapping information, quality information
		my$samples = scalar(@split) - 3 ;
		##number samples for hash
		my$sample = 0 ;
		##count sets of 3 columns for each samples
		my$count = 0 ;
		foreach my$column (3..$#split) {
			##depth
			if ($count == 0) {$count ++ ;}

			##read alignment
			elsif ($count == 1) {
				$count ++ ;
				my$depth = $split[$column-1] ;

				### Skip sites with depth less than 2 to avoid issues with binomial calculation
				if ($depth <= 1) {next ;}

				##**Call column must be processed in this order otherwise regex matches will be incorrect (e.g. all indel bp will be detected as snp calls)
				my$data = $split[$column] ;
				# Remove mapping qualities (^ marks beginning of read and the ASCII of the character following it - 33 is is the mapping quality
				$data =~ s/(\^.)//g ;
				# Remove end of read symbol
				$data =~ s/\$//g ;

				# Count the number of reads that have indels and remove indels from $data (for snp counting)
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
					#print "indel count: ", $indel, "\tindel: ", $1, "\n" ;
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
				}

				## Set absolute minimum to 5, and record minimum alternate allele count for sample+site
				if ($min_alt < 5) {$min_alt = 5 ;}
				$mpileup{$sample}{$scaffold}{$position}{"MIN_ALT"} = $min_alt ;
				#print $column, "\t", $scaffold, "\t", $position, "\t", $depth, "\t", $min_alt, "\t", $prob_accurate, "\t", $data, "\n" ;

				##record depth from previous column
				$mpileup{$sample}{$scaffold}{$position}{"DEPTH"} = $depth ;

				## Output variant call if call > min_alt at site
				$mpileup{$sample}{$scaffold}{$position}{"REF_BASE"} = $ref_base ;

				if ($ref > $min_alt) {$mpileup{$sample}{$scaffold}{$position}{"REF_COUNT"} = $ref ;}
				else {$mpileup{$sample}{$scaffold}{$position}{"REF_COUNT"} = 0 ;}

				if ($altA > $min_alt) {$mpileup{$sample}{$scaffold}{$position}{"ALT_A_COUNT"} = $altA ;}
				else {$mpileup{$sample}{$scaffold}{$position}{"ALT_A_COUNT"} = 0 ;}

				if ($altT> $min_alt) {$mpileup{$sample}{$scaffold}{$position}{"ALT_T_COUNT"} = $altT ;}
				else {$mpileup{$sample}{$scaffold}{$position}{"ALT_T_COUNT"} = 0 ;}

				if ($altC > $min_alt){$mpileup{$sample}{$scaffold}{$position}{"ALT_C_COUNT"} = $altC ;}
				else {$mpileup{$sample}{$scaffold}{$position}{"ALT_C_COUNT"} = 0 ;}

				if ($altG > $min_alt) {$mpileup{$sample}{$scaffold}{$position}{"ALT_G_COUNT"} = $altG ;}
				else {$mpileup{$sample}{$scaffold}{$position}{"ALT_G_COUNT"} = 0 ;}

				if ($indel > $min_alt) {
					$mpileup{$sample}{$scaffold}{$position}{"INDEL_COUNT"} = $indel ;
					$mpileup{$sample}{$scaffold}{$position}{"INDEL"} = $indels[0] ;
				}

				else {$mpileup{$sample}{$scaffold}{$position}{"INDEL_COUNT"} = 0 ;}

			}

			##alignment quality (not recorded)
			elsif ($count == 2) {
				##set counts for next sample
				$sample ++ ; $count = 0 ;
			}
			else {print "out of range (0-2) count: ", $count;}
		}
	}
	close PILEUP ;

	return \%mpileup ;
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
      if ( $_ =~ m/^#/ ) {next ;}
      if ( $_ =~ m/>/ ) {
        if ($seq) {$seqs{$header} = $seq ;}
        chomp ;

        my$header_line = $_ ;
				my@header_elements = split(/\s+/, $header_line) ;
        $header = $header_elements[0] ;
        $header =~ s/^>// ;
        $header =~ s/\s+$// ;

				if ($header =~ m/[|]/) {
					$header =~ s/^gi[|]\d+[|]// ;
					$header =~ s/^gb[|]// ;
					$header =~ s/[|]$// ;
					$header =~ s/^ref[|]// ;
				}

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

  foreach my$scaff (keys %seqs) {
		#delete non-keyword scaffolds
  	if ($scaff !~ m/$keyword/) {
				print $keyword, " doesn't match fasta scaffold: ", $scaff, "\n" ;
				delete $seqs{$scaff} ;
			next ;
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
