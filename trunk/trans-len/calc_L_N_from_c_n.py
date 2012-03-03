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
from nzmath.combinatorial import stirling2
from scipy import comb, factorial

def main(args):
    c, n = None, None
    try:
        opts, args = getopt.getopt(args, 'c:n:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
    for o, a in opts:
        if o == '-c':
            c = int(a)
        if o == '-n':
            n = int(a)
    if c == None or n == None:
        print >> sys.stderr, "Missing options"

if __name__ == '__main__':
    main(sys.argv[1:])

