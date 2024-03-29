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

from Bio import Entrez
import sys
import multiprocessing
import getopt
import time
import re

def print_nuc_fasta(gi):
    gb = None
    try:
        gb = Entrez.efetch(db='protein', rettype='gp', id=gi, retmode='xml')
    except Exception as err:
        print >> sys.stderr, str(err)
        pass
    print >> sys.stderr, "%s GI: %s" % (str(time.asctime()), gi)
    gb = Entrez.read(gb)[0]
    if 'GBSeq_feature-table' in gb:
        for feature_equal in gb['GBSeq_feature-table']:
            for feature in feature_equal['GBFeature_quals']:
                if feature['GBQualifier_name'] == 'coded_by':
                    region = feature['GBQualifier_value']
                    strand = 1
                    if re.match(r'complement(\S+)', region) != None:
                        strand = 2
                        region = region[11:-1]
                    region = region.split(":")
                    region[1] = region[1].split("..")
                    nuccore_id = region[0]
                    seq_start = region[1][0]
                    seq_stop = region[1][1]
                    try:
                        gi_fasta = Entrez.efetch(db='nuccore', rettype='fasta', id=nuccore_id, strand=strand, seq_stop=seq_stop, seq_start=seq_start)
                        gi_fasta = gi_fasta.read()
                        gi_fasta = gi_fasta[0] + "gi|" + gi + "|" + gi_fasta[1:]
                    except Exception as err:
                        print >> sys.stderr, str(err)
                        pass
                    return gi_fasta

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
    gis = []
    for gi_file in args:
        with open(gi_file, 'r') as fin:
            for gi in fin:
                gi = gi.strip()
                if len(gi) != 0:
                    gis.append(gi)
    pool = multiprocessing.Pool(p)
    nuc_fastas = pool.map(print_nuc_fasta, gis)
    pool.close()
    pool.join()
    for nuc_fasta in nuc_fastas:
        if nuc_fasta != None:
            print nuc_fasta            

    
    
