#!/usr/bin/python

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

import sys
from Bio import SeqIO

def CountBp(fasta_files):
    bp = 0
    for fasta_file in fasta_files:
        for seq in list(SeqIO.parse(fasta_file, 'fasta')):
            bp = bp + len(seq.seq)
    print bp
        
if __name__ == '__main__':
    CountBp(sys.argv[1:])
    sys.exit(0)

