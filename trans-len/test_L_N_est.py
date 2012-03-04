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
test my estimator for L, N using simulated data

0 < Y_{i + 1} - Y_i \leq k$
"""

import sys
import getopt

import sim_c_n
import calc_L_N_from_c_n

def main(args):
    L, N, k , runs = None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:N:k:r:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-L':
            L = int(a)
        if o == '-k':
            k = int(a)
        if o == '-N':
            N = int(a)
        if o == '-r':
            runs = int(a)
    if L == None or N == None or k == None or runs == None:
        print >> sys.stderr, "Missing options!"
        sys.exit(2)
    sim_res = sim_c_n.sim_CN(L, k, runs, N)
    LN_tab = calc_L_N_from_c_n.calc_L_N(L - 1, N, k)
    for cn_tup in sim_res['sim_CN']:
        if LN_tab[cn_tup[1]][cn_tup[0]] == None:
            print >> sys.stderr, "Error, c: %d, n: %d" % (cn_tup[1].cn_tup[0])
        else:
            # L_est, N_est
            print float(LN_tab[cn_tup[1]][cn_tup[0]][0]),
            print float(LN_tab[cn_tup[1]][cn_tup[0]][1])
    
if __name__ == '__main__':
    main(sys.argv[1:])

