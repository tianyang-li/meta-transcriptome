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
$Y_{i + 1} - Y_{i} \leq k$
"""

import getopt
import sys

import matplotlib.pyplot as plt

import help_1

def rand_read_len_prob(N, k):
    """
    cdf of contig length that can be produced
    when N reads of length k are used 
    
    it is assumed that all overlaps between reads are eqaully likely
    """
    prob = []
    for r in range(0, (N - 1) * k + 1):
        prob.append(help_1.read_split_contig(N, r, k))
    all_num = 0
    for x in prob:
        all_num += x
    prob = map(lambda x: float(x) / float(all_num), prob)
    cdf = [prob[0]]
    for i in range(len(prob) - 1):
        cdf.append(cdf[i] + prob[i + 1])
    return cdf

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
    cdf = rand_read_len_prob(N, k)
    plt.plot(range(len(cdf)), cdf)
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
