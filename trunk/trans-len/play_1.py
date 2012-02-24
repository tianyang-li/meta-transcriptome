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

import getopt
import sys
from scipy import comb

def main(args):
    n, r, k = None, None, None
    try:
        opts, args = getopt.getopt(args, 'n:r:k:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2) 
    for o, a in opts:
        if o == '-n':
            n = int(a)
        if o == '-r':
            r = int(a)
        if o == '-k':
            k = int(a)
    if n == None or r == None or k == None:
        print >> sys.stderr, "Missing arguments"
        sys.exit(2) 
    for m in range(n + 1):
        if m > r / k:
            break
        print comb(n, m) * comb(r - 1 - m * k, n - 1)
    
if __name__ == '__main__':
    main(sys.argv[1:])

