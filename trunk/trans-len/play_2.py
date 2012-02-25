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
(effective) single contig length distribution 
given (effective) length $L$, (effective) read 
length $k$, and the number of reads $N$
"""

import getopt
import sys

from nzmath.combinatorial import stirling1
from scipy import comb

import matplotlib.pyplot as plt

def single_contig_len_prb(L, N, k):
    prob = []
    #TODO
    sum = 0
    for x in prob:
        sum += x
    map(lambda x: x / sum, prob)
    return prob

def main(args):
    L, N, k = None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:N:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    if L == None or N == None or k == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    prob = single_contig_len_prb(L, N, k)
    plt.plot(range(len(prob)),prob)
    plt.grid(True)
    plt.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])


