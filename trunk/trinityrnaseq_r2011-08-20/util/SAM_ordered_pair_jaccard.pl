#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 coord_sorted.+or-.sam [windowSize=100]\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $window_length = $ARGV[1] || 100;

my $MIN_FRAGS = 100;
my $STEP_SIZE = int($window_length / 3);

main: {
	

	my @frags;

	my $sam_reader = new SAM_reader($sam_file);

	my $genome_acc = "";

	my $targeted_strand;

	my $prev_frag_end = 0;

	while ($sam_reader->has_next()) {

		my $sam_entry = $sam_reader->get_next();
		
		unless ($sam_entry->is_proper_pair()) {
			next;
		}

		unless ($sam_entry->get_mate_scaffold_name() eq "=" || $sam_entry->get_mate_scaffold_name() eq $sam_entry->get_scaffold_name()) { next; }

		unless ($sam_entry->get_aligned_position() < $sam_entry->get_mate_scaffold_position()) { next; }
		
		my $scaffold = $sam_entry->get_scaffold_name();

		unless ($genome_acc eq $scaffold) {
			## reset, switching scaffolds.
			
			if (@frags) {
				&process_frags($genome_acc, \@frags);
				@frags = ();
			}
			
			$genome_acc = $scaffold;
			$prev_frag_end = 0;
			print "variableStep chrom=$genome_acc\n";
			
		}
		
		
		my $frag_start = $sam_entry->get_aligned_position();
		my $frag_end = $sam_entry->get_mate_scaffold_position() + length($sam_entry->get_sequence());
		
		if ($prev_frag_end && $frag_start > $prev_frag_end) {
			&process_frags($genome_acc, \@frags);
			@frags = ();
			$prev_frag_end = 0;
		}

		
		push (@frags, [$frag_start, $frag_end]);
		
		#print "adding pair: $frag_start-$frag_end || $prev_frag_end\n";
		
		$prev_frag_end = $frag_end if ($frag_end > $prev_frag_end);		
	}

	if (@frags) {
		## get any remaining ones.
		&process_frags($genome_acc, \@frags);
	}

	exit(0);
}


#####
sub process_frags {
	my ($genome_acc, $frags_aref) = @_;

	my $num_frags = scalar(@$frags_aref);

	#print "// processing $num_frags frags.\n";

	unless ($num_frags > $MIN_FRAGS) { return; }
	
	my $left_range = $frags_aref->[0]->[0];
	my $right_range = $left_range;

	foreach my $coordset (@$frags_aref) {
		my ($lend, $rend) = @$coordset;
		if ($rend > $right_range) {
			$right_range = $rend;
		}
	}


	my $out_text = "";
	
	my %indices_left_window;
	my %indices_right_window;

	my $left_index = 0;
	my $right_index = 0; # track left and right indices separately.


	## compute windows
	for (my $pos = $left_range + 2*$window_length; $pos <= $right_range - 2*$window_length; $pos += $STEP_SIZE) {
		
		
		my ($left_win_a, $left_win_b) = ($pos - $window_length, $pos - int($window_length/2));
		my ($right_win_a, $right_win_b) = ($pos + int($window_length/2), $pos + $window_length);

		my $left_win = [$left_win_a, $left_win_b];
		my $right_win = [$right_win_a, $right_win_b];

		
		## empty the left window indices that are now out of range.
		foreach my $index (keys %indices_left_window) {
			if ($frags_aref->[$index]->[1] < $left_win_a) {
				delete $indices_left_window{$index};
				
			}
		}

		## empty the right window indices that are now out of range
		foreach my $index (keys %indices_right_window) {
			if ($frags_aref->[$index]->[1] < $right_win_a) {
				delete $indices_right_window{$index};
			}
		}

		
		## walk through frags till encounter one in range
		while ($left_index <= $#$frags_aref 
			   && 
			   $frags_aref->[$left_index]->[0] <= $left_win_b) {
			
			if (&overlaps($left_win, $frags_aref->[$left_index])) {
				$indices_left_window{$left_index} = 1;
			}

			$left_index++;
		}
		
		## 
		# walk through right frags from last pointer.
		while ($right_index <= $#$frags_aref
			   &&
			   $frags_aref->[$right_index]->[0] <= $right_win_b) {
			
			if (&overlaps($right_win, $frags_aref->[$right_index])) {
				$indices_right_window{$right_index} = 1;
			}

			$right_index++;
		}
		
		if (%indices_left_window && %indices_right_window) {
			my $jaccard = &compute_jaccard(\%indices_left_window, \%indices_right_window);
			
			$out_text .= "$pos\t" . sprintf("%.3f", $jaccard) . "\n";
			
		}
	}
	
	if ($out_text) {
		print $out_text;
	}
	
	return;

}


####
sub overlaps {
	my ($coordset_a, $coordset_b) = @_;

	if ($coordset_a->[0] < $coordset_b->[1]
		&&
		$coordset_a->[1] > $coordset_b->[0]) {
		return(1); # yes, overlaps
	}
	else {
		return(0);
	}
}

####
sub compute_jaccard {
	my ($left_index_href, $right_index_href) = @_;

	my %any;
	my $both_count = 0;

	foreach my $index (keys %$left_index_href, keys %$right_index_href) {
		$any{$index} = 1;
	}

	my @indices = keys %any;
	my $sum_count = scalar(@indices);

	foreach my $index (@indices) {
		if ($left_index_href->{$index} && $right_index_href->{$index}) {
			$both_count++;
		}
	}

	

	my $jaccard_coeff = $both_count/$sum_count;
	return($jaccard_coeff);
}

