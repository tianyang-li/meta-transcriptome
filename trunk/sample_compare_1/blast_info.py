#!/usr/bin/python

from Bio.Blast import NCBIXML

import sys
import itertools

blastrec = NCBIXML.parse(open(sys.argv[1]))

for rec in blastrec:
    if len(rec.alignments) > 0:
        for (desc, al) in itertools.izip(rec.descriptions, rec.alignments):
            raw_input("################## Press a key to continue...............")
