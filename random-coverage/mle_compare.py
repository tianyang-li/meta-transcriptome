#!/usr/bin/env python

import sys
import getopt

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
    return est_L

def mle_L(d, N):
    est_L = []
    for l in range(d + 1):
        est = int(float(N * l) / float(N - 1))
        if est <= l:
            est = l + 1
        else:
            if (est - l) * ((est + 1) ** N) <= (est + 1 - l) * (est ** N):
                est += 1
        est_L.append(est)
    return est_L

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
    umve = umve_L(d, N)
    mle = mle_L(d, N)
    for Lu, Lm in zip(umve, mle):
        print Lu, Lm, Lu - Lm
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
