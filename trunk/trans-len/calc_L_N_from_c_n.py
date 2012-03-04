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
import json

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
        if n==N:
            sum1 = L - c
        else:
            return 0
    return sum1 * comb(N, n, exact=True) * help_1.read_split_contig(n, c, k)

def calc_L_N(c, n, k):
    """
    @return: a 2D list containing estimators for all $(c, n)$, 
             entry [c, n] is the corresponding [L, N]
    """
    LN_tab = [[None, [Fraction(1), Fraction(1)]]]
    if c == 0:
        for n_val in range(2, n + 1):
            L, N = 1, n_val
            tot_cn_num = 0
            L_tmp1, N_tmp1 = 0, 0
            for n_tmp in range(1, N):
                tmp_cn_num = get_cn_num(L, N, 0, n_tmp, k)
                tot_cn_num += tmp_cn_num
                L_tmp1 += (tmp_cn_num * LN_tab[0][n_tmp][0])
                N_tmp1 += (tmp_cn_num * LN_tab[0][n_tmp][1])
            tmp_cn_num = get_cn_num(L, N, 0, N, k)
            tot_cn_num += tmp_cn_num
            L_est = (tot_cn_num - L_tmp1) / tmp_cn_num
            N_est = (tot_cn_num * N - N_tmp1) / tmp_cn_num
            LN_tab[0].append([L_est, N_est])
    else:
        for n_val in range(1, n + 1):
            L, N = c + 1, n_val
            tot_cn_num = 0
            L_tmp1, N_tmp1 = 0, 0
            for c_prev in range(c):
                for n_prev in range(1, N + 1):
                    tmp_cn_num = get_cn_num(L, N, c_prev, n_prev, k)
                    tot_cn_num += tmp_cn_num
                    L_tmp1 += (tmp_cn_num * LN_tab[c_prev][n_prev][0])
                    N_tmp1 += (tmp_cn_num * LN_tab[c_prev][n_prev][1])
            for n_prev in range(1, N):
                tmp_cn_num = get_cn_num(L, N, c, n_prev, k)
                tot_cn_num += tmp_cn_num
                L_tmp1 += (tmp_cn_num * LN_tab[c][n_prev][0])
                N_tmp1 += (tmp_cn_num * LN_tab[c][n_prev][1])
            tmp_cn_num = get_cn_num(L, N, c, N, k)
            tot_cn_num += tmp_cn_num
            L_est = (tot_cn_num * L - L_tmp1) / tmp_cn_num
            N_est = (tot_cn_num * N - N_tmp1) / tmp_cn_num
            if N == 1:
                LN_tab.append([None, [L_est, N_est]])
            else:
                LN_tab[c].append([L_est, N_est])
    return LN_tab

def main(args):
    c, n, k = None, None, None
    try:
        opts, args = getopt.getopt(args, 'c:n:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-c':
            c = int(a)
        if o == '-n':
            n = int(a)
        if o == '-k':
            k = int(a)
    if c == None or n == None or k == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    LN_tab = calc_L_N(c, n, k)
    float_LN_tab = []
    for c_list in LN_tab:
        for LN_tup in c_list:
            if LN_tup == None:
                float_LN_tab.append([None])
            else:
                float_LN_tab[-1].append([float(LN_tup[0]), float(LN_tup[1])])
    print json.dumps(float_LN_tab)
    

if __name__ == '__main__':
    main(sys.argv[1:])

