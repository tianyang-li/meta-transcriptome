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

import sys
import random
import multiprocessing
import getopt
import json

class SimSingleParams(object):
    def __init__(self, L, k, N, runs):
        self.L = L
        self.k = k
        self.N = N
        self.runs = runs

def sim_single(params):
    L = params.L  #effective length
    k = params.k
    N = params.N
    runs = params.runs
    single_contig_lens = []
    random.seed()
    for i in range(runs):
        read_pos = []
        for j in range(N):
            read_pos.append(random.randint(0, L - 1))       
        single_contig_len = SingleContigLen(read_pos, k)
        if single_contig_len != 0:
            single_contig_lens.append(single_contig_len)
    return single_contig_lens

def SingleContigLen(read_pos, k):
    read_pos.sort()
    contig_len = 0
    pos = read_pos[0]
    for next_pos in read_pos[1:]:
        if next_pos != pos:
            if pos + k <= next_pos:
                return 0
            contig_len += (next_pos - pos)
            pos = next_pos
    contig_len += k
    return contig_len

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'L:N:k:r:p:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    L = None  
    k = None
    N = None
    runs = None
    p = 1  # number of parallel computations
    for o, a in opts:
        if o == '-L':
            L = int(a)
        elif o == '-k':
            k = int(a)
        elif o == '-r':
            runs = int(a)
        elif o == '-N':
            N = int(a)
        elif o == '-p':
            p = int(a)
    if (L == None) or (k == None) or (N == None) or (runs == None):
        print >> sys.stderr, "Missing options!"
        sys.exit(2)
    L = L - k + 1  #effective length
    chunksize = runs / p + 1
    
    pool = multiprocessing.Pool(processes=p)
    sim_single_params = []
    for i in range(p):
        sim_single_params.append(SimSingleParams(L, k, N, chunksize))
    results = pool.map(sim_single, sim_single_params)
    pool.close()
    pool.join()
    contig_lens = []
    for result in results:
        contig_lens.extend(result)
    L = L + k - 1
    len_distr = {}
    runs = chunksize * p
    for i in range(1, L + 1):
        len_distr[i] = float(contig_lens.count(i)) / float(runs)
    sim_result = {"L": L, "k": k, "N": N, "contig_len": len_distr}
    print json.dumps(sim_result)
    

