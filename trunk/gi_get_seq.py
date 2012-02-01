#!/usr/bin/env python

# Copyright (C) 2011 Tianyang Li
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
./gi-get-seq.py file-containing-GI output-seq-file
"""

from Bio import Entrez
import sys
import multiprocessing
import getopt

def print_nuc_fasta(params):
    gi = params[0]
    lock = params[1]
    gb = Entrez.efetch(db='protein', rettype='gb', id=gi, retmode='xml')
    gb = Entrez.read(gb)
    with lock:
        # TODO: GenBank get nuc fasta 
        print gb

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'p:m:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    p = 1
    email = "example@example.com"
    for o, a in opts:
        if o == '-p':
            p = int(a)
        elif o == '-m':
            email = a
    Entrez.email = email
    mgr = multiprocessing.Manager()
    gis = []
    for gi_file in args:
        with open(gi_file, 'r') as fin:
            for gi in fin:
                gi = gi.strip()
                if len(gi) != 0:
                    gis.append([gi, mgr.Lock()])
    chunksize = len(gis) / p
    pool = multiprocessing.Pool(p)
    pool.map(print_nuc_fasta, gis, chunksize=chunksize)
    pool.close()
    pool.join()
                
    
    
