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

"""
test uniformity of reads distribution on contigs using the 
Rao's spacing test
"""

import sys
import getopt
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from Bio import SeqIO
import HTSeq
import json

def main(args):
    bowtie2, read_len, contigs, fout = None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'b:c:l:o:')
    except getopt.GetoptError as err:
        print >> sys.stderr(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-b':
            # bowtie2 SAM file
            bowtie2 = arg
        if opt == '-l':
            read_len = int(arg)
        if opt == '-c':
            # assembled contigs (FASTA)
            contigs = arg
        if opt == '-o':
            # output file
            fout = arg
    if bowtie2 == None or read_len == None or contigs == None or fout == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])



