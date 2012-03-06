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
from HTSeq import SAM_Reader
import json
from scipy.constants import pi

def main(args):
    bowtie2, read_len, contigs, fout, alpha = None, None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'b:c:l:o:a:')
    except getopt.GetoptError as err:
        print >> sys.stderr(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            # alpha used in test
            alpha = float(arg)
        if opt == '-b':
            # bowtie2 SAM file
            bowtie2 = arg
        if opt == '-l':
            read_len = int(arg)
        if opt == '-c':
            # assembled contigs (FASTA)
            contigs = arg
        if opt == '-o':
            # output file containing ids that do not pass the test
            fout = arg
    if bowtie2 == None or read_len == None or contigs == None or fout == None or alpha == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
    
    contigs_len = {}
    for rec in SeqIO.parse(contigs, 'fasta'):
        contigs_len[rec.id] = [len(rec.seq), []]

    for align in SAM_Reader(bowtie2):
        if align.aligned:
            contigs_len[align.iv.chrom][1].append(align.iv.start)

    with open(fout, 'w') as res:
        circular = importr('circular')
        rst = circular.rao_spacing_test
        circ = circular.circular
        
        for contig_id, contig_entry in contigs_len.items(): 
            if contig_entry[0] >= read_len:
                eff_len = contig_entry[0] - read_len + 1
                pos_in_rad = []
                for start_pos in contig_entry[1]:
                    if start_pos < eff_len:
                        pos_in_rad.append(2 * pi * float(start_pos) / float(eff_len))
                if len(pos_in_rad) > 3:
                    rst_res = rst(circ(robj.FloatVector(pos_in_rad)), alpha=alpha)
                    rst_res = str(rst_res)
                    print rst_res
                    rst_res = rst_res.split("\n")[5]
                    if "Reject" in rst_res:
                        # id, length, # of reads
                        res.write("%s %d %d\n" % (contig_id, contig_entry[0], len(contig_entry[1])))
                
    
if __name__ == '__main__':
    main(sys.argv[1:])



