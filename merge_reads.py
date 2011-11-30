#!/usr/bin/python

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
./merge_reads.py [fasta/fastq] [output] [input files]
"""

import sys
from Bio import SeqIO

def main(argv):
    reads = []
    for fin in argv[3:]:
        for read in SeqIO.parse(fin, argv[2]):
            reads.append(read)
    SeqIO.write(reads, argv[2], argv[1])

if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)

