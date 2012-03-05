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
from Bio import SeqIO

def main(args):
    db, b6 = None, None
    try:
        opts, args = getopt.getopt(args, 'd:b:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-d':
            # file containing all the db seqs, assumed to be in FASTA
            db = a
        if o == '-b':
            # usearch b6 files saying which reads mapped to db seqs
            b6 = a
    if db == None or b6 == None:
        print >> sys.stderr, "missing arguments"
        sys.exit(2)

    db_ids = set([])
    for rec in SeqIO.parse(db, 'fasta'):
        db_ids.add(rec.id)
        
    with open(b6,'r') as b6_fin:
        for entry in b6_fin:
            entry = entry.strip().split("\t")
            contig_id = entry[1].split(" ")[0]
            read_id = entry[0]
            if contig_id in db_ids:
                print read_id


if __name__ == '__main__':
    main(sys.argv[1:])

