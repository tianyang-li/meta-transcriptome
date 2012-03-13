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

import sys
import getopt
from Bio import SeqIO
from HTSeq import SAM_Reader

def complete_cover(cover, read_len):
    if cover[0] == 0:
        return False 
    if len(cover) - read_len >= 0:
        if cover[len(cover) - read_len] == 0:
            return False
    prev = None
    for pos in xrange(len(cover)):
        if cover[pos] != 0:
            if prev != None:
                if pos - prev >= read_len:
                    return False
            prev = pos
    return True

def main(args):
    contigs, sam, read_len = None, None, None
    try:
        opts, args = getopt.getopt(args, 'c:s:l:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
    for opt, arg in opts:
        if opt == '-c':
            contigs = arg
        if opt == '-s':
            sam = arg
        if opt == '-l':
            read_len = int(arg)
    if contigs == None or sam == None:
        print >> sys.stderr, "Missing options"
    
    contigs_cover = {}
    for rec in SeqIO.parse(contigs, 'fasta'):
        if len(rec.seq) >= read_len:
            contigs_cover[rec.id] = [0] * len(rec.seq)
    
    for align in SAM_Reader(sam):
        if align.aligned:
            if align.iv.chrom in contigs_cover:
                contigs_cover[align.iv.chrom][align.iv.start] += 1
    
    for contig_id, contig_cover in contigs_cover.iteritems():
        if complete_cover(contig_cover, read_len):
            print contig_id
        
    
if __name__ == '__main__':
    main(sys.argv[1:])



