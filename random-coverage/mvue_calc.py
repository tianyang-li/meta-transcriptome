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

def umve_L(d, N):
    """
    d -- X_max - X_min
    N -- number of samples
    @return: estimated L (umve)
    """
    est_L = []
    for l in range(d + 1):
        est = (l + 1) ** (N + 1)
        for x, prev_L in zip(range(len(est_L)), est_L):
            if x == 0:
                est -= (l + 1)
            else:
                est -= prev_L * (l + 1 - x) * ((x + 1) ** N + (x - 1) ** N - 2 * (x ** N))
        if l != 0:
            est = float(est) / float((l + 1) ** N - 2 * (l ** N) + (l - 1) ** N)
        est_L.append(est)
    return est_L[-1]

def main(args):
    N = None
    d = None
    try:
        opts, args = getopt.getopt(args, 'd:N:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-d':
            d = int(a)
        if o == '-N':
            N = int(a)
    if (N == None) or (d == None):
        print >> sys.stderr, "Missing options!"
        sys.exit(2)
    est_L = umve_L(d, N)
    print est_L
    
if __name__ == '__main__':
    main(sys.argv[1:])
