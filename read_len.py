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
./read_len.py [fasta/fastq] [read files]
"""

import sys
from Bio import SeqIO

def main(argv):
    for fin in argv[2:]:
        for read in SeqIO.parse(fin, argv[1]):
            print read.id, len(str(read.seq))

if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)

