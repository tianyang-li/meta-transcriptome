#!/usr/bin/env python

import sys
import getopt

from help_1 import read_split_contig as split

def main(args):
    n, c, d = None, None, None
    try:
        opts, args = getopt.getopt(args, 'n:c:d:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-n':
            n = int(arg)
        if opt == '-c':
            c = int(arg)
        if opt == '-d':
            d = int(arg)
    if n == None or c == None or d == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    
    for n_ in xrange(1, n):
        for c_ in xrange(c):
            good_val = split(n_, c_, d)
            test_val = (c_ + 1) * split(n_ - 1, c_, d)
            print test_val
            for i in xrange(1, min(d, c_) + 1):
                test_val += split(n_ - 1, c_ - i, d)
                print n_ - 1, c_ - i, test_val
            print n_, c_, good_val, test_val
            if test_val != good_val:
                print >> sys.stderr, "you are wrong"
                sys.exit(1)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])    


