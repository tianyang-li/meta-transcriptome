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
remove reads that were aligned >= 2 times and contigs associated with them
"""

import sys
import getopt
from HTSeq import SAM_Reader
from Bio import SeqIO

def main(args):
    contigs, sam, fout_prefix, reads, fmt = None, None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'c:s:o:r:f:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-c':
            contigs = arg
        if opt == '-s':
            sam = arg
        if opt == '-o':
            fout_prefix = arg
        if opt == '-r':
            reads = arg
        if opt == '-f':
            fmt = arg
    if contigs == None or sam == None or fout_prefix == None or reads == None:
        print >> sys.stderr, "missing options"
        sys.exit(1)
        
    read_aln_count = {}
    for rec in SeqIO.parse(reads, fmt):
        read_aln_count[rec.id] = []    
    contig_aln_ok = {}
    for rec in SeqIO.parse(contigs, 'fasta'):
        contig_aln_ok[rec.id] = True
    
    for aln in SAM_Reader(sam):
        if aln.aligned:
            read_aln_count[aln.read.name].append(aln.iv.chrom)
    
if __name__ == '__main__':
    main(sys.argv[1:])    



