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

import sys
import getopt
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from Bio import SeqIO
import datetime

import contig_bowtie2_uniform_test

def main(args):
    read_len, seq_db = None, None
    try:
        opts, b6s = getopt.getopt(args, 'l:s:')
    except  getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-l':
            read_len = int(a)
        if o == '-s':
            seq_db = a
    if read_len == None or seq_db == None or b6s == []:
        print >> sys.stderr, "Missing options or arguments!"
        sys.exit(2)
    db_len = contig_bowtie2_uniform_test.get_fasta_len(seq_db)
    

if __name__ == '__main__':
    main(sys.argv[1:])


