#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 jaccard.wig\n\n";

my $jaccard_file = $ARGV[0] or die $usage;

my $WIN_VALLEY_SEARCH = 100;
my $WIN_HILL_SEARCH = 300;

my $VALLEY_MAX = 0.05;
my $HILL_MIN = 0.25;


my $VERBOSE = 0;

main: {
	
	my %scaff_to_positions = &parse_wig_entries($jaccard_file);
		
	foreach my $scaff (sort keys %scaff_to_positions) {
		
		print "variableStep chrom=$scaff\n";

		my @vals = @{$scaff_to_positions{$scaff}};
		
		for (my $i = 0; $i <= $#vals; $i++) {
			
			my ($coord, $val) = @{$vals[$i]};
			
			#print "-visiting: $coord, $val\n";

			if ($val <= $VALLEY_MAX) {
				
				print "Exploring valley: $coord, $val\n" if $VERBOSE;
				
				my ($min_pos, $min_val) = &explore_valley($i, $coord, $val, \@vals);
				
				
				my ($coord, $val) = @{$vals[$min_pos]};
				
				print "\tgot valley low at: $coord, $val\n" if $VERBOSE;
				
				if (&neighboring_hill($min_pos, \@vals)) {
					
					
					print "$coord\t1\n";

					my $placeholder = $coord;
					while ($coord < $placeholder + $WIN_VALLEY_SEARCH 
						   && 
						   $i <= $#vals) {
						
						$i++;
						
						my $aref = $vals[$i];
						unless (ref $aref) {
							print STDERR "warning, no aref for $i [ignore for now]\n";
							last;
						}
						
						($coord, $val) = @{$vals[$i]};
						print "Skipping through: $coord, $val\n" if $VERBOSE;
					}
				}
				

			}
		}
	}

	
	
	exit(0);
}

#### 
sub parse_wig_entries {
	my ($jaccard_file) = @_;

	my %scaffold_to_entries;

	my $curr_scaff;
	
	open (my $fh, $jaccard_file) or die "Error, cannot open file $jaccard_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		chomp;
		if (/^variableStep chrom=(\S+)/) {
			$curr_scaff = $1;
			next;
		}
			
		my ($coord, $val) = split(/\t/);
		
		unless (defined($coord) && defined($val)) {
			#die "Error, line $_ lacks expected format";
			next;
		}
		
		push (@{$scaffold_to_entries{$curr_scaff}}, [$coord, $val]);
	}

	close $fh;


	return(%scaffold_to_entries);
}

####
sub explore_valley {
	my ($pos, $valley_coord, $min_val, $vals_aref) = @_;

	
	my $smallest_val = $min_val;
	my $smallest_val_pos = $pos;
	
	## look right:
	for (my $i = $pos+1; $i <= $#$vals_aref; $i++) {
		my ($coord, $val) = @{$vals_aref->[$i]};
		if ($coord > $valley_coord + $WIN_VALLEY_SEARCH) {
			last;
		}
		if ($val < $smallest_val) {
			$smallest_val = $val;
			$smallest_val_pos = $i;
		}
	}
	
	
	return($smallest_val_pos, $smallest_val);

}

####
sub neighboring_hill {
	my ($pos, $vals_aref) = @_;

	print "\t\texploring neighboring hill from $pos\n" if $VERBOSE;
	
	#print "$pos\t$vals_aref\n";
	
	my ($start_coord, $start_val) = @{$vals_aref->[$pos]};

	## search left
	my $found_hill_left_flag = 0;
	
	for (my $i = $pos-1; $i >= 0; $i--) {
		my ($coord, $val) = @{$vals_aref->[$i]};
		print "\t\t\t\t$coord, $val\n" if $VERBOSE;

		if ($coord < $start_coord - $WIN_HILL_SEARCH) {
			last;
		}

		if ($val >= $HILL_MIN) {
			$found_hill_left_flag = 1;
			last;
		}
	}

	if ($found_hill_left_flag) {
		print "\t\t\t\tFound Left Hill.\n" if $VERBOSE;
	}

	unless ($found_hill_left_flag) {
		return(0); 
	}

	## search right
	for (my $i = $pos + 1; $i <= $#$vals_aref; $i++) {
		my ($coord, $val) = @{$vals_aref->[$i]};
		print "\t\t\t\t$coord, $val\n" if $VERBOSE;
		if ($coord > $start_coord + $WIN_HILL_SEARCH) {
			last;
		}

		if ($val >= $HILL_MIN) {
			print "\t\t\t\tFound Right Hill.\n" if $VERBOSE;
			return(1); # got it!
		}
	}

	print "\t\t\t\tNo right hill\n" if $VERBOSE;
	
	return(0); # no hill
}


