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
    for i in range(L):
        sim_cn_tab.append([0] * (N + 1))
    for cn_tup in sim_res['sim_CN']:
        sim_cn_tab[cn_tup[1]][cn_tup[0]] += 1
    sim_cn_tab = map(lambda c_list: map(lambda cnt: float(cnt) / float(len(sim_res['sim_CN'])), c_list), sim_cn_tab)
    print "####"
    for c_list in sim_cn_tab:
        for prob in c_list:
            print prob,
        print ""
    
    calc_cn_tab = []
    calc_sum = 0
    for c in range(L):
        calc_cn_tab.append([0])
        for n in range(1, N + 1):
            tmp_cnt = calc_L_N_from_c_n.get_cn_num(L, N, c, n, k)
            calc_cn_tab[c].append(tmp_cnt)
            calc_sum += tmp_cnt
    calc_cn_tab = map(lambda c_list: map(lambda cnt: float(cnt) / float(calc_sum), c_list), calc_cn_tab)
    print "####" 
    for c_list in calc_cn_tab:
        for prob in c_list:
            print prob,
        print ""   

if __name__ == '__main__':
    main(sys.argv[1:])

