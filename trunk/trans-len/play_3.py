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
given $N$ reads and (effective) read length $k$
"""

import getopt
import sys

from nzmath.combinatorial import stirling1
from scipy import comb

import matplotlib.pyplot as plt

def rand_read_len_prob(N, k):
    prob = []
    #TODO
    sum = 0
    for x in prob:
        sum += prob
    map(lambda x: x / sum, prob)
    return prob

def main(args):
    N, k = None, None
    try:
        opts, args = getopt.getopt(args, 'N:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-N':
            N = int(a)
        if o == '-k':
            k = int(a)
    if N == None or k == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    prob = rand_read_len_prob(N, k)
    plt.plot(range(len(prob)), prob)
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
