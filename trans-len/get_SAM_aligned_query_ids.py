#!/usr/bin/env python

#  Copyright (C) 2012 Tianyang Li
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License

"""
get query ids for queries algined to db seqs with db ids
"""

import sys
import getopt
from HTSeq import SAM_Reader

def main(args):
    sam, db_ids = None, None
    try:
        opts, args = getopt.getopt(args, 's:d:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-s':
            sam = arg
        if opt == '-d':
            db_ids = arg
    if sam == None or db_ids == None:
        print >> sys.stderr, "missing arguments"
        sys.exit(2)
    
    dids = set([])
    with open(db_ids, 'r') as fin:
        for line in fin:
            dids.add(line.strip())
    
    for align in SAM_Reader(sam):
        if align.aligned:
            if align.iv.chrom in dids:
                print align.read.name

if __name__ == '__main__':
    main(sys.argv[1:])


