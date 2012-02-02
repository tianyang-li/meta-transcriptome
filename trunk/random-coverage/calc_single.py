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
    def __init__(self, L, k, n):
        self.L = L
        self.k = k
        self.n = n
        self.contig_len = contig_len

def calc_single(params):
    L = params.L
    k = params.k
    n = params.n
    prob_list = []
    #TODO: calculate this probability
    for cl in range(1, L + 1):
        prob = 0
        prob_list.append(prob)
    return prob_list

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
    calc_params = []
    for n in range(1, L + 1):  # n is the number of possible starting positions
        calc_params.append(CalcSingleParams(L, k, n))
    pool = multiprocessing.Pool(p)
    prob_lists = pool.map(calc_single, calc_params)
    pool.close()
    pool.join()
    probs = {}
    for n in range(1, L + 1):
        probs[n] = 0
    for n, prob_list in zip(range(1, L + 1), prob_lists):
        for cl, prob in zip(range(1, L + 1), prob_list):
            probs[cl] += prob
    calc_result = {'L': L, 'k': k, 'N': N, 'contig_len': probs}
    print json.dumps(calc_result)
    
        
