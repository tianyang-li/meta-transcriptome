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
filter out "random" contigs 
"""

import sys
import getopt
from Bio import SeqIO
from math import exp

def keep_contig(c, n, d):
    """
    c - effective contig length (# of possible read starting positions)
    n - number of reads on the contig
    d - maximum distance between consecutive read starting positions
    
    @return: True if contig doesn't seem random 
    """
    if n <= 2:
        return False
    # expression level: lambda (Poisson parameter)
    pparam = n / (c + 2 * d)
    # prob of no reads
    p0 = exp(-pparam)
    # prob when there are reads
    p1 = 1 - p0
    # expected length
    el = 0.0
    for i in range(1, d + 1):
        el += (p1 * i)
        p1 *= p0
    if ((n - 1) * el > c):
        return True
    else:
        return False
    

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
    
    contigs_len = {}
    for rec in SeqIO.parse(contigs, 'fasta'):
        if len(rec.seq) >= read_len:
            contigs_len[rec.id] = [len(rec.seq), 0]
    
    for b6 in b6s:
        with open(b6, 'r') as fin:
            for entry in fin:
                entry = entry.strip().split("\t")
                contig_id = entry[1].split(" ")[0]
                contigs_len[contig_id][1] += 1
    
    for contig_id, contig_entry in contigs_len.items():
        if keep_contig(contig_entry[0] - read_len + 1, contig_entry[1], read_len - kmer + 1):
            print contig_id

if __name__ == '__main__':
    main(sys.argv[1:])




