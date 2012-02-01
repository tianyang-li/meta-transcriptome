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

import sys
import scipy
import multiprocessing
import getopt

def main(sysargs):
    try:
        opts, args = getopt.getopt(sysargs, 'L:N:k:r:p:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    L = None # effective length
    k = None
    N = None
    runs = None
    p = 1
    for o, a in opts:
        if o == '-L':
            eL = int(a)
        elif o == '-k':
            k = int(a)
        elif o == '-r':
            runs = int(a)
        elif o == '-N':
            N = int(a)
        elif o == '-p':
            p = int(a)
    if (eL == None) or (k == None) or (N == None) or (runs == None):
        print >> sys.stderr, "Missing options!"
        sys.exit(2)
    L = L - k + 1 
    
if __name__ == '__main__':
    main(sys.argv[1:])
