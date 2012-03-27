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

from __future__ import division

import sys
import getopt

import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def llss(c, n, L, N, d_max):
    """
    log likelihood (sufficient statistic)
    !!! it is assumed that the likelihood is never 0 !!!
    """
    if n == N:
        return math.log(L - c) - N * math.log(L)
    ll = 0.0
    for i in xrange(N - n, N + 1):
        ll += math.log(i)
    ll -= (N * math.log(L))
    pos_cnt = 0
    cur_avail_pos = L - 2 * d_max - c - 1
    if cur_avail_pos > 0:
        pos_cnt = (L - 2 * d_max - c) * (cur_avail_pos ** (N - n))
    for i in xrange(d_max):
        cur_avail_pos = L - i - d_max - c - 1
        if cur_avail_pos > 0:
            pos_cnt += (2 * (cur_avail_pos ** (N - n)))
    ll += math.log(pos_cnt)
    return ll


def main(args):
    c, n , d_max = None, None, None
    try:
        opts, args = getopt.getopt(args, 'c:n:d:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-c':
            c = int(arg)
        if opt == '-n':
            n = int(arg)
        if opt == '-d':
            d_max = int(arg)
    if c == None or n == None or d_max == None:
        print >> sys.stderr, "missing options"
        sys.exit(1)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    Ns, Ls, lls = [], [], []
    
    # TODO: calculate likelihood
    
    Ns = np.array(Ns)
    Ls = np.array(Ls)
    lls = np.array(lls)
    
    ax.scatter(Ns, Ls, lls)
    
    ax.set_xlabel('N')
    ax.set_ylabel('L')
    ax.set_zlabel('log likelihood')
    
    plt.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])    




