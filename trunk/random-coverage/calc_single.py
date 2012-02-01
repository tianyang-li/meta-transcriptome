#!/usr/bin/env python

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

import getopt
import sys
import multiprocessing
import json
import scipy

class CalcSingleParams(object):
    def __init__(self, L, k, N, contig_len):
        self.L = L
        self.k = k
        self.N = N
        self.contig_len = contig_len

def calc_single(params):
    L = params.L
    k = params.k
    N = params.N
    cl = params.contig_len
    prob = 0
    #TODO: calculate this probability
    return cl, prob

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'L:N:k:p:')
    except getopt.GetoptError as err:
            print >> sys.stderr, str(err)
            sys.exit(2)        
    p = 1
    L = None
    k = None
    N = None
    for o, a in opts:
        if o == '-L':
            L = int(a)
        elif o == '-N':
            N = int(a)
        elif o == '-k':
            k = int(a)
    if (N == None) or (L == None) or (k == None):
        print sys.stderr, "Missing options"
        sys.exit(2)
    chunksize = (L - k + 1) / p 
    calc_params = []
    for contig_len in range(k, L + 1):
        calc_params.append(CalcSingleParams(L, k, N, contig_len))
    pool = multiprocessing.Pool(p)
    results = pool.map(calc_single, calc_params, chunksize=chunksize)
    pool.close()
    pool.join()
        
