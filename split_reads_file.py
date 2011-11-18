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

"""
split_reads_file.py [input type] [input file] 
    [# of reads in each splitted file] [output prefix]
"""

from Bio import SeqIO
import sys

def main(argv):
    each_reads = []
    k = 0
    split_count = 0
    for read in SeqIO.parse(argv[2], argv[1]):
        each_reads.append(read)
        if k < int(argv[3]):
            k += 1
        else:
            k = 0
            SeqIO.write(each_reads, argv[4] + (".%d.%s" % (split_count, argv[1])), argv[1])
            split_count += 1
            each_reads = []
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)

