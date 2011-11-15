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
#
#  You should have received a copy of the GNU General Public License

# TODO: add paired end read analysis

"""
Analysis metatranscriptome sequences using de Bruijn graph
"""

__authors__ = [
               '"Tianyang Li" <tmy1018@gmail.com>',
]

import getopt
import sys

import dbg_metatrans

def Usage(prog):
    print >> sys.stderr, "\nUsage:"
    print >> sys.stderr, "%s [-k kmer-length] [-s single-read-FASTQ1] [-s single-read-FASTQ2] ..." % prog
    print >> sys.stderr, "\nOptions:\n"
    print >> sys.stderr, "    -h    Print this short help message"
    print >> sys.stderr, "    -k    Kmer length"
    print >> sys.stderr, "    -s    A FASTQ files containing single reads"
    print >> sys.stderr, ""
    
def main(argv):
    single = []  # single read FASTQ    
    opts, args = getopt.getopt(argv[1:], "hk:s:")
    required_opts = 0
    for o, a in opts:
        if o == "-h":
            Usage(argv[0])
            sys.exit(1)
        elif o == "-k":
            k = int(a)
            required_opts += 1         
        elif o == "-s":
            single.append(a)
            required_opts += 1
        else:
            print >> sys.stderr, "Unrecognized option: %s" % o
            print >> sys.stderr, "See\n    %s -h\nfor more information" % argv[0]
            sys.exit(1)
    if required_opts < 2:
        print >> sys.stderr, "Too few options and arguments"
        print >> sys.stderr, "See\n    %s -h\nfor more information" % argv[0]       
    dbg_metatrans.dbg_metatrans(k, single)

if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)
