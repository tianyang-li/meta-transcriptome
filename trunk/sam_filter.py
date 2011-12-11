#!/usr/bin/env python

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
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Author Tianyang Li
#  E-Mail tmy1018 gmail com

"""
filter reads according to alignment in SAM file
Usage:
./sam_filter.py [0/1] [SAM] [fasta/fastq] [read] [output]
0 - keep those not aligned in SAM
1 - keep those aligned in SAM
"""

import sys
import HTSeq
from Bio import SeqIO

def main(argv):
    keep = True
    if argv[1] == '0':
        keep = False
    if argv[1] == '1':
        keep = True
        
    aligned = set([])
        
    for align in HTSeq.SAM_Reader(argv[2]):
        if align.aligned:
            aligned.add(align.read.name)
    
    filtered = []
    for read in SeqIO.parse(argv[4], argv[3]):
        if keep:
            if read.name in aligned:
                filtered.append(read)
        else:
            if read.name not in aligned:
                filtered.append(read)
    
    SeqIO.write(filtered, argv[5], argv[3])

if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)

