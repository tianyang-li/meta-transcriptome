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

'''
calculate nucleotide frequencies in FASTA files
'''

import sys
from Bio import SeqIO

def FastaNFreq(fasta_files):
    bp = 0
    a = 0
    t = 0
    c = 0
    g = 0
    for fasta_file in fasta_files:
        for seq in list(SeqIO.parse(fasta_file, 'fasta')):
            bp = bp + len(seq.seq)
            for nuc in seq.seq:
                if nuc == 'A' or nuc == 'a':
                    a = a + 1
                elif nuc == 'C' or nuc == 'c':
                    c = c + 1
                elif nuc == 'G' or nuc == 'g':
                    g = g + 1
                elif nuc == 'T' or nuc == 't':
                    t = t + 1
    print "A:\t%f" % (float(a) / float(bp))
    print "T:\t%f" % (float(t) / float(bp))
    print "C:\t%f" % (float(c) / float(bp))
    print "G:\t%f" % (float(g) / float(bp))

if __name__ == '__main__':
    FastaNFreq(sys.argv[1:])
    sys.exit(0)
