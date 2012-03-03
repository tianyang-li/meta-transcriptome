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
compute estimators for all $(contig_len, contig_reads)$
$0 \leq contig_len \leq c$, $1 \leq contig_reads \leq n$
"""

import getopt
import sys
from scipy import comb
from fractions import Fraction

import help_1

def get_cn_num(L, N, c, n, k):
    sum1 = 0
    if L - (2 * k + c + 1) >= 0:
        if L - (2 * k + c + 1) == 0:
            sum1 += 1
        else:
            sum1 += (((L - (2 * k + c + 1)) ** (N - n)) * (L - (2 * k + c)))
    for i in range(k):
        if L - (i + k + c + 1) >= 0:
            if L - (i + k + c + 1) == 0:
                sum1 += 2
            else:
                sum1 += (2 * (L - (i + k + c)) * ((L - (i + k + c + 1)) ** (N - n)))
    if sum1 == 0:
        return 0
    return sum1 * comb(N, n, exact=True) * help_1.read_split_contig(n, c, k)

def main(args):
    c, n, k = None, None, None
    try:
        opts, args = getopt.getopt(args, 'c:n:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
    for o, a in opts:
        if o == '-c':
            c = int(a)
        if o == '-n':
            n = int(a)
        if o == '-k':
            k = int(a)
    if c == None or n == None or k == None:
        print >> sys.stderr, "Missing options"

if __name__ == '__main__':
    main(sys.argv[1:])

