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
exact multinomial test

-b bowtie2 results
-l read length (assumed to be all equal)
-c assembled contigs
"""

import sys
import getopt
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from Bio import SeqIO
import HTSeq
import datetime

def count_bowtie2_align(bowtie2, db_len):
    for align in HTSeq.SAM_Reader(bowtie2):
        if align.aligned:
            if db_len[align.iv.chrom][0] > align.iv.start and align.iv.start >= 0:
                db_len[align.iv.chrom][1][align.iv.start] += 1
                db_len[align.iv.chrom][2] += 1
        
        
def get_fasta_len(fasta_file):
    """
    get length of each entry in a fasta file,
    return it as a dict 
    """
    fasta_len = {}
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        # [length of seq (int), reads distribution (lsit), total reads count]
        fasta_len[str(rec.id)] = [len(rec.seq), [0] * len(rec.seq), 0]
    return fasta_len

def main(args):
    bowtie2, read_len, contigs = None, None, None
    try:
        opts, args = getopt.getopt(args, 'b:l:c:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-b':
            bowtie2 = a
        if o == '-l':
            read_len = int(a)
        if o == '-c':
            contigs = a
    if bowtie2 == None or read_len == None or contigs == None:
        print >> sys.stderr, "Missing options!"
        sys.exit(2)
    contigs_len = get_fasta_len(contigs)
    count_bowtie2_align(bowtie2, contigs_len)
    EMT = importr('EMT')
    mt = EMT.multinomial_test
    for contig_entry, i in zip(contigs_len.values(), range(1, len(contigs_len) + 1)):
        if contig_entry[2] > 2:
            print >> sys.stderr, "#%d: length: %d, reads: %d, time: %s" % (i, contig_entry[0], contig_entry[2], datetime.datetime.now())
            pv = mt(robj.IntVector(contig_entry[1][:-(read_len - 1)]), robj.FloatVector([1 / float(contig_entry[0] - read_len + 1)] * (contig_entry[0] - read_len + 1)), MonteCarlo=True, ntrial=1000000)            
            # contig length, p-value
            print contig_entry[0], pv[-1][0]
    
if __name__ == '__main__':
    main(sys.argv[1:])
