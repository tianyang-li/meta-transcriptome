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

import sys
import getopt

import calc_L_N_from_c_n
import sim_c_n
from numpy import array
from scipy.stats import chisquare, chi2

def main(args):
    runs, L, N, k = None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'r:L:N:k:')
    except getopt.GetoptError as err:
        print sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-r':
            # # of runs
            runs = int(a)
        if o == '-L':
            L = int(a)
        if o == '-N':
            N = int(a)
        if o == '-k':
            k = int(a)
    if runs == None or L == None or N == None or k == None:
        print sys.stderr, "Missing options"
        sys.exit(2)
    
    sim_res = sim_c_n.sim_CN(L, k, runs, N)
    sim_cn_tab = []
    sim_cn_tot = 0
    for i in range(L):
        sim_cn_tab.append([0] * (N + 1))
    for cn_tup in sim_res['sim_CN']:
        sim_cn_tab[cn_tup[1]][cn_tup[0]] += 1
        sim_cn_tot += 1
    
    calc_cn_tab = []
    calc_sum = 0
    for c in range(L):
        calc_cn_tab.append([0])
        for n in range(1, N + 1):
            tmp_cnt = calc_L_N_from_c_n.get_cn_num(L, N, c, n, k)
            calc_cn_tab[c].append(tmp_cnt)
            calc_sum += tmp_cnt
    calc_cn_tab = map(lambda c_list: map(lambda cnt: float(cnt) / float(calc_sum), c_list), calc_cn_tab)
    
    obs_freq = []
    exp_freq = []
    for sim_c, calc_c in zip(sim_cn_tab, calc_cn_tab):
        for sim_freq, calc_prob in zip(sim_c, calc_c):
            if calc_prob != 0:
                obs_freq.append(sim_freq)
                exp_freq.append(calc_prob * sim_cn_tot)
    print obs_freq
    print exp_freq
    chi_stat, pv = chisquare(array(obs_freq), f_exp=array(exp_freq), ddof=(len(exp_freq) - 1))
    pv = 1 - chi2.cdf(chi_stat, len(exp_freq) - 1)
    print pv

if __name__ == '__main__':
    main(sys.argv[1:])

