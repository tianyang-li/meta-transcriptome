#!/usr/bin/env python

# Copyright (C) 2012 Tianyang Li
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

"""
test the MLE for estimating contig length and total number of reads
"""

import sys
import getopt

def main(args):
    L, N, runs, d_max = None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:N:r:d:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-L':
            # effective contig length
            L = int(arg)
        if opt == '-N':
            # total number of reads 
            N = int(arg)
        if opt == '-d':
            # max difference between 2 different read starting positions
            d_max = int(arg)
        if opt == '-r':
            # number of runs when generating random contigs
            runs = int(arg)
    if L == None or N == None or runs == None or d_max == None:
        print >> sys.stderr, "missing options"
        sys.exit(1)
    
if __name__ == '__main__':
    main(sys.argv[1:])    



