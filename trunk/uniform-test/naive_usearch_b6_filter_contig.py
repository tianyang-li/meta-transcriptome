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
    read_len, kmer, contigs = None, None, None
    try:
        opts, b6s = getopt.getopt(args, 'r:k:c:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-c':
            # name of the file that contain assembled contigs
            contigs = a
        if o == '-k':
            # kmer length used in assembly
            kmer = int(a)
        if o == '-r':
            # read length
            read_len = int(a)
    if contigs == None or kmer == None or read_len == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    

if __name__ == '__main__':
    main(sys.argv[1:])




