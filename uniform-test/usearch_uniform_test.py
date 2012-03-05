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
import json

import contig_bowtie2_uniform_test

def count_b6_map(db_len, b6s):
    for b6 in b6s:
        with open(b6, 'r') as fin:
            for entry in fin:
                entry = entry.strip().split("\t")
                db_seq = entry[1].split(" ")[0]
                start_pos = int(entry[8])
                db_len[db_seq][2] += 1
                db_len[db_seq][1][start_pos - 1] += 1
                

def main(args):
    read_len, seq_db, ntrial, fout = None, None, None, None
    try:
        opts, b6s = getopt.getopt(args, 'l:s:t:o:')
    except  getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-l':
            # read length, assumed to be all equal
            read_len = int(a)
        if o == '-s':
            # sequence database
            seq_db = a
        if o == '-t':
            # # of trials used in monte carlo multinomial test
            ntrial = int(a)
        if o == '-o':
            # output file
            fout = a
    if read_len == None or seq_db == None or b6s == [] or ntrial == None or fout == None:
        print >> sys.stderr, "Missing options or arguments!"
        sys.exit(2)
        
    db_len = contig_bowtie2_uniform_test.get_fasta_len(seq_db)
    count_b6_map(db_len, b6s)
    
    fout = open(fout, 'w')
    EMT = importr('EMT')
    mt = EMT.multinomial_test
    
    for seq_entry in db_len.values():
        if seq_entry[2] > 2 and seq_entry[0] >= read_len:
            print >> sys.stderr, "####################"
            print >> sys.stderr, datetime.datetime.now(), " len: ", seq_entry[0], " reads: ", seq_entry[2]
            print >> sys.stderr, "####################"
            # test uniformity on the whole contig
            pv_all = mt(robj.IntVector(seq_entry[1]), robj.FloatVector([1 / float(seq_entry[0])] * seq_entry[0]), MonteCarlo=True, ntrial=ntrial)
            # test uniformity on the effective part of the contig
            pv_eff = mt(robj.IntVector(seq_entry[1][:-(read_len - 1)]), robj.FloatVector([1 / float(seq_entry[0] - read_len + 1)] * (seq_entry[0] - read_len + 1)), MonteCarlo=True, ntrial=ntrial)
            uni_test = {'len': seq_entry[0], 'reads': seq_entry[2], 'read_pos': seq_entry[1]}
            uni_test['all_test'] = pv_all[-1][0]
            uni_test['eff_test'] = pv_eff[-1][0]
            fout.write("%s\n" % json.dumps(uni_test))
            
    fout.close()

if __name__ == '__main__':
    main(sys.argv[1:])


