#!/usr/bin/env perl

# Copyright (C) 2011 Tianyang Li
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

# simluate uniform random covering of a
# transcript by reads of equal length
#
# this is a model for Illumina single end reads

use strict;
use warnings;

use Getopt::Std;

# -L	length of transcript
# -N	number of reads
# -k 	length of read
# -t	number of threads
my %params = ( t => "1" );
getopts( 'L:N:k:t:', \%params );

# calculate probability using my formula
{
	use bignum;

	sub calc_prob {
		my $params = shift;
		my $k      = $params->{k};
		my $N      = $params->{N};
		my $L      = $params->{L} - $k + 1;    # effective length
		my $prob = 0;
		$prob;
	}

	printf "%g\n", &calc_prob( \%params );
}

