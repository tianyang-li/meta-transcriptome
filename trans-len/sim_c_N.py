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
-r    number of runs for each possible $N$
-N    maximum $N$ to use in simulation (choose this according to Lander-Waterman?)
"""

import sys
import getopt
import random
import os
import base64
import matplotlib.pyplot as plt

def main(args):
    L, k, runs, N_max = None, None, None, None
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
            N_max = int(a)
    if L == None or k == None or runs == None or N_max == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    seed = os.urandom(64)
    random.seed(seed)
    seed_b64 = base64.b64encode(seed)
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])


