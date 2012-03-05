#!/usr/bin/env python

import sys
import getopt

from calc_L_N_from_c_n import get_cn_num as cnt
import sim_c_n

def bf_cnt_cn(L, N, k):
    cn_cnt = []
    for i in range(L ** N):
        start_pos = [0] * L
        for e in range(N):
            start_pos[i % L] += 1
            i /= L
        cn_cnt.extend(sim_c_n.count_CN(start_pos, k))
    return cn_cnt

def main(args):
    L, N, k = None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:N:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-L':
            L = int(a)
        if o == '-N':
            N = int(a)
        if o == '-k':
            k = int(a)
    if L == None or N == None or k == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    bf_cn = bf_cnt_cn(L, N, k)
    cn_cnt = []
    for c in range(L):
        cn_cnt.append([[0]] * (N+1))
    for c in range(L):
        for n in range(1,N+1):
            cn_cnt[c][n].append(cnt(L, N, c, n, k))
    for cn_tup in bf_cn:
        cn_cnt[cn_tup[1]][cn_tup[0]][0] += 1
    for c_list in cn_cnt:
        for x in c_list:
            print x,
        print ""
            

if __name__ == '__main__':
    main(sys.argv[1:])
