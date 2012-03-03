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
simulate tuples of $(c, N, \pm 1)$ where
$c$ is the (effective) contig length and 
$N$ is the number of reads on the contig
and the result is $1$ if this contig is not the
only one from the transcript, $-1$ if the contig is 
the only one from the transcript

-L    (effective) transcript length
-k    max distance between 2 read starting positions, (effective) read length
          $Y_{i + 1} - Y_i \leq k$
-r    number of runs 
-N    number of reads on transcript 
"""

import sys
import getopt
import random
import json
import os
import base64

def count_CN(start_pos, k):
    cn_tups = []
    n_c, l_c = 0, 0
    prev = None
    for pos in range(len(start_pos)):
        if start_pos[pos] != 0:
            if prev == None:
                n_c = start_pos[pos]
            else:
                if pos - prev <= k:
                    n_c += start_pos[pos]
                    l_c += pos - prev
                else:
                    cn_tups.append([n_c, l_c, 1])
                    n_c = start_pos[pos]
                    l_c = 0
            prev = pos
    if n_c != 0:
        cn_tups.append([n_c, l_c, 1])
    if len(cn_tups) == 1:
        cn_tups[0][2] = -1
    return cn_tups

def sim_CN(L, k, runs, N):
    """
    simulate tuples of $(c, N, \pm 1)$ where
    $c$ is the (effective) contig length and 
    $N$ is the number of reads on the contig
    and the result is $1$ if this contig is not the
    only one from the transcript, $-1$ if the contig is 
    the only one from the transcript
    
    L: (effective) transcript length
    k: max distance between 2 read starting positions, (effective) read length
        $Y_{i + 1} - Y_i \leq k$
    runs: number of runs
    N: number of reads on transcript
    
    @return: a dictionary containing base64 encoded seed and tuples
    """
    seed = os.urandom(64)
    random.seed(seed)
    seed_b64 = base64.b64encode(seed)
    
    cn_tups = []
    for r in range(runs):
        start_pos = [0] * L
        for i in range(N):
            start_pos[random.randint(0, L - 1)] += 1
        cn_tups.extend(count_CN(start_pos, k))
    
    return {'seed_b64': seed_b64, 'sim_CN': cn_tups}

def main(args):
    L, k, runs, N = None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:k:r:N:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-L':
            L = int(a)
        if o == '-k':
            k = int(a)
        if o == '-r':
            runs = int(a)
        if o == '-N':
            N = int(a)
    if L == None or k == None or runs == None or N == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    sim_res = sim_CN(L, k, runs, N)
    for res in sim_res['sim_CN']:
        for x in res:
            print x,
        print ""
    

if __name__ == '__main__':
    main(sys.argv[1:])


