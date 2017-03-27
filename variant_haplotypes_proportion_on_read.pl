use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;
use Statistics::Distributions qw< chisqrprob > ;

#Usage: perl test_reads_for_haplotypes.pl pval positions2test.txt sample.bam

### NOTES: 
######### Optical duplicate reads are filtered out
######### No need to look for reverse complement matches - they are already complemented to the reference and the flag is marked 16 for reverse complement
######### The remaining duplicate reports for read matches at a site are real - some reads overlap (e.g. PJRI47 ILLUMINA-D00365:314:HA5YEADXX:2:1116:2462:4203) - These are accounted for at line 119-142: if 2 reads in pair report the same site, they are counted as separate reads

# Chi-squared test determines if variant sites present on the same read represent a significant association,
# or could have been produced by change (e.g. if allele frequencies are high for both sites).
# Haplotype associations are inferred by the presence of more variants on one read/pair than expected by chance.
# Expected values for x variants on a given read were calculated as the product of the probabilities of selecting each independently
# E = V1_count/depth * V2_count/depth * V3_count/depth ...
# Observed values are calculated by finding reads where more than one variant site is present and reporting how many of these occurred on the read
# O = variants_on_read/variants_in_read_range 
# 
# Degrees of freedom = #variant sites on read - 1 
# Dof justification: If "the number of degrees of freedom is the number of values in the final calculation of a statistic that are free to vary.",
# then since there are x sites that can only be in one other state, off the read, there are x values, and for chi-squared tests, dof = x values - 1 

my$threshold = $ARGV[0] ;
my$positions = $ARGV[1] ;
my$bam = $ARGV[2] ;

# Output: read	scaffold	#variants_in_range #variants_on_read	proportion	expected_value	chisq	pval
my$output1 = basename($positions) ;
$output1 =~ s/AF.txt/read_co-occurring_variants.txt/ ;
$output1 =~ s/AF_counts.txt/read_co-occurring_variants.txt/ ;
$output1 =~ s/all_AF_counts_subpopsegsites.txt/subpopsegsites_read_co-occurring_variants.txt/ ;
$output1 =~ s/.filtered// ;

# Output: scaffold	haplotype	#reads_supporting 
my$output2 = basename($positions) ;
$output2 =~ s/AF.txt/reads_supporting_haplotypes.txt/ ;
$output2 =~ s/AF_counts.txt/reads_supporting_haplotypes.txt/ ;
$output1 =~ s/all_AF_counts_subpopsegsites.txt/subpopsegsites_reads_supporting_haplotypes.txt/ ;
$output2 =~ s/.filtered// ;

open POSITIONS, "<$positions" ;

my%positions ;

# Get positions to test
# Only using minor count because only considering biallelic, so only have 2 haplotype states (associated variants (2 haplotypes) vs. not)
while (<POSITIONS>) {
	chomp ;
	my@split = split(/\t/, $_) ;	
	my$scaffold = $split[0] ;
	my$pos = $split[1] ;
	my$ref_base = $split[2] ;
	my$ref_count = $split[3] ;
	my$alt_base = $split[4] ;
	my$alt_count = $split[5] ;
	my$allele_freq = $split[6] ;
	
	if ($alt_count < $ref_count) {
		$positions{$scaffold}{$pos}{"NT"} = $alt_base ;	
		$positions{$scaffold}{$pos}{"MINOR_COUNT"} = $alt_count ;
	}
	
	if ($ref_count < $alt_count) {
		$positions{$scaffold}{$pos}{"NT"} = $ref_base ;
		$positions{$scaffold}{$pos}{"MINOR_COUNT"} = $ref_count ;
	}

	$positions{$scaffold}{$pos}{"DEPTH"} = $ref_count + $alt_count ;
}

close POSITIONS ;

open BAM, "samtools view $bam |" ;

my%reads ;
my%haplotypes ;
my$dups = 0 ;
my$processed_reads = 0 ;

# Get reads that map to those positions
while (<BAM>) {
	next if (/^(\@)/) ;
	s/\n// ;
	s/\r// ;
	my@sam = split(/\t/, $_) ;

	my$read = $sam[0] ;
	my$flag = $sam[1] ;
	my$scaffold = $sam[2] ;

	if ($scaffold eq "Sv_mito_chromosome") {next ;}
	# Skip optical/PCR duplicates
	my@dup_flags = (1024, 1040, 1089, 1097, 1105, 1107, 1113, 1123, 1137, 1153, 1161, 1171, 1177, 1185, 1187, 1201) ;
	if (grep {$_ == $flag} @dup_flags) {
		$dups ++ ;
		next;
	}
	
	my$results ;
	if ($sam[5] =~ m/[ID]/) {$results = read_with_indel($_, \%positions) ;}
	else {$results = read_no_indel($_, \%positions) ;}

	my%results = %{$results} ;	

	my$variants_read_range = $results{"READ_RANGE"} ;
	my$variants_on_read = $results{"ON_READ"} ;
	my@expected = @{$results{"EXPECTED"}} ;
	my@haplotype= @{$results{"HAPLOTYPES"}} ;

	my$string = "" ;

	## If read contains a variant site, record it in the hash (or append if read pair was recorded)
	if ($variants_read_range >= 1) {
		if (exists $reads{$read}{$scaffold}) {
		
			$reads{$read}{$scaffold}{"variants_read_range"} += $variants_read_range ;
			$reads{$read}{$scaffold}{"variants_on_read"} += $variants_on_read ;

			# Some reads overlap, making the sites reported more than once. 
			# Eliminate extra accounts of haplotype from calculation below ($haplotype[$i] aligned to $expected[$i])
			# Subtract a variant_read_range and variant_on_read for each duplicate haplotype call removed
			# Only appending novel sites:
			foreach my$i (0..$#haplotype) {
				if (grep {$_ != $haplotype[$i]} @{$reads{$read}{$scaffold}{"haplotype_sites"}}) {
					$reads{$read}{$scaffold}{"expected_value"} *= $expected[$i] ;
					push @{$reads{$read}{$scaffold}{"haplotype_sites"}}, $haplotype[$i] ;
					$string .= "$haplotype[$i]\t" ;
				}
				
				else {
					$reads{$read}{$scaffold}{"variants_read_range"} -- ;
					$reads{$read}{$scaffold}{"variants_on_read"} -- ;
				}	
			}

		}
		
		else {
			$reads{$read}{$scaffold}{"variants_read_range"} = $variants_read_range ;
			$reads{$read}{$scaffold}{"variants_on_read"} = $variants_on_read ;
			$reads{$read}{$scaffold}{"haplotype_sites"} = \@haplotype ;
			$reads{$read}{$scaffold}{"expected_value"} = 1 ;
			foreach my$i (0..$#expected) {
				$reads{$read}{$scaffold}{"expected_value"} *= $expected[$i] ;
				$string .= "$haplotype[$i]\t" ;
			}
		
		}
		
		if (exists $haplotypes{$scaffold}{$string}) {$haplotypes{$scaffold}{$string} ++ ;}
		else {$haplotypes{$scaffold}{$string} = 1 ;}
		

	}
	
	$processed_reads ++ ;
	
}	

close BAM ;

open OUT1, ">$output1" or die "cannot open $output1" ;
open OUT2, ">$output2" or die "cannot open $output2" ;

my$total = 0 ;
my$all_on_off_read = 0 ;
my$significant = 0 ;

foreach my$read (nsort keys %reads) {

	my@scaffs = keys %{$reads{$read}} ;
	if (scalar @scaffs > 1) {
		print $read, "\n" ;
	}
	my@haplotypes0 = @{$reads{$read}{$scaffs[0]}{"haplotype_sites"}} ;
	
	if (scalar @scaffs == 2) {
		my$variants_in_range = $reads{$read}{$scaffs[0]}{"variants_read_range"} + $reads{$read}{$scaffs[1]}{"variants_read_range"} ;
		my$variants_present = $reads{$read}{$scaffs[0]}{"variants_on_read"} + $reads{$read}{$scaffs[1]}{"variants_on_read"} ;
		my$proportion = $variants_present/$variants_in_range ;
		
		if ($variants_present < 2) {next;}
		
		my$expected = $reads{$read}{$scaffs[0]}{"expected_value"} * $reads{$read}{$scaffs[1]}{"expected_value"} ;
		my@haplotypes1 = @{$reads{$read}{$scaffs[1]}{"haplotype_sites"}} ;
		
		#chi-squared test, df = sites-1
		my$chisq = ($proportion - $expected)*($proportion - $expected) / $expected ;
		my$df = $variants_in_range - 1 ;
		my$pval = chisqrprob($df, $chisq) ;
		
		print OUT1 $read, "\t", $scaffs[0], "-", $scaffs[1], "\t", join(",", @haplotypes0), "-", join(",", @haplotypes1), "\t", $variants_in_range, "\t", $variants_present, "\t", $proportion, "\t", $expected, "\t", $chisq, "\t", $pval, "\n" ;
		
		if ($pval <= $threshold) {$significant ++ ;}
		
		if ($proportion == 1 || $proportion == 0) {$all_on_off_read ++ ;}
	}
	
	else {				
		if ($reads{$read}{$scaffs[0]}{"variants_on_read"} < 2) {next;}
		
		my$proportion = $reads{$read}{$scaffs[0]}{"variants_on_read"}/$reads{$read}{$scaffs[0]}{"variants_read_range"} ;
		
		#chi-squared test, df = sites-1
		my$chisq = ($proportion - $reads{$read}{$scaffs[0]}{"expected_value"})*($proportion - $reads{$read}{$scaffs[0]}{"expected_value"}) / $reads{$read}{$scaffs[0]}{"expected_value"} ;
		my$df = $reads{$read}{$scaffs[0]}{"variants_read_range"} - 1 ;
		my$pval = chisqrprob($df, $chisq) ;

		print OUT1 $read, "\t", $scaffs[0], "\t", join(",", @haplotypes0), "\t", $reads{$read}{$scaffs[0]}{"variants_read_range"}, "\t", $reads{$read}{$scaffs[0]}{"variants_on_read"}, "\t", $proportion, "\t", $reads{$read}{$scaffs[0]}{"expected_value"}, "\t", $chisq, "\t", $pval, "\n" ;

		if ($pval <= $threshold) {$significant ++ ;}
		
		if ($proportion == 1 || $proportion == 0) {$all_on_off_read ++ ;}
		
	}		

	$total ++ ;

}

my%haplotype_counts ;

foreach my$scaffold (nsort keys %haplotypes) {
	foreach my$haplotype (nsort keys %{$haplotypes{$scaffold}}) {
		my@haplotypes = split(/\t/, $haplotype) ;
		
		if (scalar(@haplotypes) > 1) {
			print OUT2 $scaffold, "\t", $haplotype, "\t", $haplotypes{$scaffold}{$haplotype}, "\n" ;
			
			foreach (@haplotypes) {$haplotype_counts{$_} = () ;}
		}
	}
}

close OUT1 ;
close OUT2 ;

my$haplotype_counts = scalar keys %haplotype_counts ;

print "Fraction of reads with complete haplotype:  ", $all_on_off_read, " / ", $total, " = ", $all_on_off_read/$total, "\n" ;
print "Fraction of reads with significantly associated haplotype:  ", $significant, " / ", $total, " = ", $significant/$total, "\n" ;
print "Sites in a linked haplotype: ", $haplotype_counts, "\n" ;
print "Total reads processed: ", $processed_reads, "\n" ;
print "Optical/PCR duplicate reads removed: ", 	$dups, "\n" ;


sub read_no_indel {
	my@sam = split(/\t/, $_[0]) ;
	my%positions = %{$_[1]} ;
		
	my$read = $sam[0] ;
	my$scaffold = $sam[2] ;
	my$read_start = $sam[3] ;
	my$sequence = $sam[9] ;	
	
	# Search read sequence for variant sites
	my@sequence = split("", $sequence) ;
	my$read_end = $read_start + length($sequence)-1 ;
	my@range = $read_start .. $read_end ;
	
	my$variants_read_range = 0 ;
	my$variants_on_read = 0 ;
	my@haplotype = () ;
	
	my@expected = () ;
	
	foreach my$i (0..$#range) {
		if ($positions{$scaffold}{$range[$i]}{"NT"}) {
			$variants_read_range ++ ;
				
			if ($sequence[$i] eq $positions{$scaffold}{$range[$i]}{"NT"}) {
				push @expected, $positions{$scaffold}{$range[$i]}{"MINOR_COUNT"} / $positions{$scaffold}{$range[$i]}{"DEPTH"} ;
				$variants_on_read ++ ;
				push @haplotype, $range[$i] ;
			}
		}
	}

	my%results ;
	$results{"READ_RANGE"} = $variants_read_range ;
	$results{"ON_READ"} = $variants_on_read ;
	$results{"EXPECTED"} = \@expected ;
	$results{"HAPLOTYPES"} = \@haplotype ;
	
	my%haplotypes = map{$_, 1} @haplotype ;
	if (keys %haplotypes < scalar(@haplotype)) {
		print "keys: ", scalar keys %haplotypes, "\thaplotypes: ", scalar(@haplotype), "\tnt: ", join(",", @haplotype), "\n" ;
	} 

	return \%results ;	
}

sub read_with_indel {
	my@sam = split(/\t/, $_[0]) ;
	my%positions = %{$_[1]} ;

	my$read = $sam[0] ;
	my$flag = $sam[1] ;
	my$scaffold = $sam[2] ;
	my$read_start = $sam[3] ;
	my$cigar = $sam[5] ;
	my$sequence = $sam[9] ;

	## Parse indels (add/subtract them from postions to recreate reference sequence length)
	my@cigars = () ;
	my@insertions = () ;
	my@deletions = () ;
	if ($cigar =~ m/[ID]/) {	
		while ($cigar =~ m/([0-9]+[MIDS=X])/g) {push @cigars, $1 ;}	
		
		my$read_counter = 0 ;
		my$new_sequence ;
		foreach my$cigar (@cigars) {
			my$bp = $cigar ;
			$bp =~ s/([MIDS=X])//g ;
			my$type = $1 ;
			if ($type eq "M" || $type eq "=" || $type eq "S" || $type eq "X") {
				$new_sequence .= substr($sequence, $read_counter, $bp) ;
				$read_counter += $bp ;
			}

			elsif ($type eq "I") {
				#record indel
				push @insertions, substr($sequence, $read_counter, $bp) ;

				#skip adding indel sequence to new sequence to generate reference length 
				$read_counter += $bp ;				
			}
			
			else {
				#add Ds to new sequence for deleted region
				foreach my$i (1..$bp) {
					$new_sequence .= "D" ;
				}
				push @deletions, $bp ;
				#don't increment counter because deletions aren't in read sequence ;
			}
		}
	
		$sequence = $new_sequence ;
		
	}
	
	# Search read sequence for variant sites	
	my@sequence = split("", $sequence) ;
	my$read_end = $read_start + length($sequence)-1 ;
	my@range = $read_start .. $read_end ;
	
	my$variants_read_range = 0 ;
	my$variants_on_read = 0 ;
	my@haplotype = () ;
	
	my@expected = () ;
	
	my$cigar_counter = 0 ;

## What is the variant type? Insertion, deletion, or snp?
## Approach: iterate through sequence bp by bp checking the cigar string to know how to deal with variation
## Matches are checked for snps, insertions are checked for identical indel sequence, and deletions are checked for start site and length relative to the reference
	foreach my$i (0..$#range) {
		#iterate through cigar string too
#		if ($cigar =~ m/[ID]/) {
#			print "Read position: ", $i+1, "\t", "Cigar position: ", $cigar_counter, "\tCounter gap: ", $i+1 - $cigar_counter, "\t", "Cigar string: ", join("", @cigars), "\t", "Cigar index 0: ", $cigars[0], "\n" ;
#		}
				
		if ($cigars[0] && $cigar_counter < $i+1) {

			if ($cigars[0] =~ m/[IMS=X]/) {
				## Check matching sequence first
				if ($cigars[0] =~ m/[MS=X]/) {
					## check for snp call
					if ($positions{$scaffold}{$range[$i]}{"NT"}) {
						my$variant = $positions{$scaffold}{$range[$i]}{"NT"} ;
						$variants_read_range ++ ;

						if ($variant && $sequence[$i] eq $variant) {
#							print $read, "\n", "snp: ", $scaffold, ":", $range[$i], "\t", $variant, "\t", $sequence[$i], "\n", $sam[9], "\n", $sequence, "\n\n" ;
							push @expected, $positions{$scaffold}{$range[$i]}{"MINOR_COUNT"} / $positions{$scaffold}{$range[$i]}{"DEPTH"} ;
							$variants_on_read ++ ;
							push @haplotype, $range[$i] ;
						}
					}
				
					## check for end of match to delete cigar
					my$bp = $cigars[0] ;
					$bp =~ s/[MIDS=X]// ;
					if (($i+1 - $cigar_counter) == $bp) {
#						print "cigar bp: ", $bp, "\t", "Deleted: ", $cigars[0], "\n" ;
						$cigar_counter += $bp ;
						@cigars = splice(@cigars, 1) ;
					}
				}
				
				## Immediately take care of insertions to remove them from the cigar string
				## Unless indel is at the start of the sequence, you have to subset by $i+1 because this is technically reading the cigar string from one position ahead, which has been cut from the read sequence
				if ($cigars[0] && $cigars[0] =~ m/I/) {
					## check for insertion call
					if ($positions{$scaffold}{$range[$i]}{"NT"}) {
						my$variant = $positions{$scaffold}{$range[$i]}{"NT"} ;
						$variants_read_range ++ ;
						
						my$read_insert ;
						if ($variant && $cigar =~ m/^$cigars[0]/) {$read_insert = substr($sam[9], $i, length($variant)-1) ;}
						else { if ($variant) { $read_insert = substr($sam[9], $i+1, length($variant)-1) ;}}
						
						if ($read_insert eq substr($variant,1)) {
							push @expected, $positions{$scaffold}{$range[$i]}{"MINOR_COUNT"} / $positions{$scaffold}{$range[$i]}{"DEPTH"} ;
							$variants_on_read ++ ;
#							print "INSERTION", "\t", $scaffold, "\t", $range[$i], "\n", $read, "\n", $sam[9], "\n", $sequence, "\n\n" ;
							push @haplotype, $range[$i] ;
						}
					}
						
					## delete insertion cigar
					my$bp = $cigars[0] ;
					$bp =~ s/[MIDS=X]// ;
#					print "cigar bp: ", $bp, "\t", "Deleted: ", $cigars[0], "\n" ;
					@cigars = splice(@cigars, 1) ;
				}
			}
			
			elsif ($cigars[0] =~ m/D/) {
				## check for deletion call
				if ($positions{$scaffold}{$range[$i]}{"NT"}) {
					my$variant = $positions{$scaffold}{$range[$i]}{"NT"} ;
					$variants_read_range ++ ;

					my$bp = $cigars[0] ;
					$bp =~ s/D// ;

					# deletion is a match if its position and length are the same as in the variant call (valid b/c relative to the reference, this info is unambiguous)
					if ($variant && $cigar_counter == $i+1 && $bp == length($variant)-1) {
						push @expected, $positions{$scaffold}{$range[$i]}{"MINOR_COUNT"} / $positions{$scaffold}{$range[$i]}{"DEPTH"} ;
						$variants_on_read ++ ;
#						print "DELETION", "\t", $scaffold, "\t", $range[$i], "\n", $read, "\n", $sam[9], "\n", $sequence, "\n\n" ;
						push @haplotype, $range[$i] ;
					}
				}	
				
				## delete deletion cigar
				my$bp = $cigars[0] ;
				$bp =~ s/[MIDS=X]// ;

				if (($i+1 - $cigar_counter) == $bp) {
#					print "cigar bp: ", $bp, "\t", "Deleted: ", $cigars[0], "\n" ;
					$cigar_counter += $bp ;
					@cigars = splice(@cigars, 1) ;
				}
			}
			
			else {
				print "Cigar $cigars[0] not matched\n" ;
			}
		}
	}
	
	my%results ;
	$results{"READ_RANGE"} = $variants_read_range ;
	$results{"ON_READ"} = $variants_on_read ;
	$results{"EXPECTED"} = \@expected ;
	$results{"HAPLOTYPES"} = \@haplotype ;

	return \%results ;	
	
	if (keys %haplotypes < scalar(@haplotype)) {
		print "keys: ", scalar keys %haplotypes, "\thaplotypes: ", scalar(@haplotype), "\tnt: ", join(",", @haplotype), "\n" ;
	} 

}