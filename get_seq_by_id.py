#!/usr/bin/python

#  Copyright (C) 2011 Tianyang Li
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
following the conventions of Biopython,
get sequence records from files using rec.id
"""

import getopt
import sys
from Bio import SeqIO

def main(args):
    ids, fmt = None, None
    try:
        opts, seq_fins = getopt.getopt(args, 'i:f:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-i':
            ids = a
        if o == '-f':
            fmt = a
    if ids == None or fmt == None or seq_fins == []:
        print >> sys.stderr, "missing arguments"
        sys.exit(2)
        
    seq_ids = set([])
    with open(ids, 'r') as ids_fin:
        for entry in ids_fin:
            entry = entry.strip()
            seq_ids.add(entry)
    
    for seq_fin in seq_fins:
        for seq_rec in SeqIO.parse(seq_fin, fmt):
            if seq_rec.id in seq_ids:
                print seq_rec.format(fmt)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])


