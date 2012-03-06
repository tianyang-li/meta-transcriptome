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
from Bio import SeqIO

def main(args):
    fmt = None
    try:
        opts, args = getopt.getopt(args, 'f:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-f':
            fmt = arg
    if fmt == None:
        print >> sys.stderr, "Missing format"
        sys.exit(2)
    
    seq_cnt = 0
    for fin in args:
        for rec in SeqIO.parse(fin, fmt):
            seq_cnt += 1
    print seq_cnt


if __name__ == '__main__':
    main(sys.argv[1:])
    
    
