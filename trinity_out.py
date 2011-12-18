#!/usr/bin/env python

#  Copyright (C) 2011 Tianyang Li
#
#  tmy1018 (at) gmail (dot) com
#
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
Usage:
./trinity_out.py [Trinity.fasta] [read SAM]
"""

import sys
import string
from Bio import SeqIO
import HTSeq

def main(argv):
    contigs = {}
    for frag in SeqIO.parse(argv[1], 'fasta'):
        contigs[frag.name] = [len(frag.seq), 0]
    for align in HTSeq.SAM_Reader(argv[2]):
        if align.aligned:
            contigs[align.iv.chrom][1] += len(align.read)
    for contig in contigs.keys():
        print float(contigs[contig][1]) / float(contigs[contig][0])

if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)


