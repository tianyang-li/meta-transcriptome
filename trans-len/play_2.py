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
$Y_{i + 1} - Y_{i} \leq k$
"""

import getopt
import sys

import matplotlib.pyplot as plt

import help_1

def single_contig_len_prb(L, N, k):
    prob = []
    for r in range(L):
        prob.append((L - r) * help_1.read_split_contig(N, r, k))
    all_num = 0
    for x in prob:
        all_num += x
    prob = map(lambda x: float(x) / float(all_num), prob)
    cdf = [prob[0]]
    for i in range(len(prob) - 1):
        cdf.append(cdf[i] + prob[i + 1])
    return prob

def main(args):
    L, N, k = None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:N:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-L':
            L = int(a)
        if o == '-N':
            N = int(a)
        if o == '-k':
            k = int(a)
    if L == None or N == None or k == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    cdf = single_contig_len_prb(L, N, k)
    plt.plot(range(len(cdf)), cdf)
    plt.grid(True)
    plt.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])


