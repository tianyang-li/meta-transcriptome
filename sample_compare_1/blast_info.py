#!/usr/bin/python

from Bio.Blast import NCBIXML

import sys
import itertools

blastrec = NCBIXML.parse(open(sys.argv[1]))

for rec in blastrec:
    if len(rec.alignments) > 0:
        print "*************************************** Query START"
        print rec.query
        print "*************************************** Query END"
        for (desc, al) in itertools.izip(rec.descriptions, rec.alignments):
            print "Description"
            print "title: %s" % desc.title
            print "e: %lf" % desc.e
            print "number of alignments: %d" % desc.num_alignments
            print "Alignment"
            print "title: %s" % al.title
            print "HSPs:"
            for hsp in al.hsps:
                print "\t%s" % hsp.query
                print "\t%s" % hsp.match
                print "\t%s" % hsp.sbjct
            raw_input("################## Press a key to continue...............")
