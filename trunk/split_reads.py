#!/usr/bin/env python

from Bio import SeqIO
import sys
import getopt

def main(args):
    try:
        opts, args = getopt.getopt(args, 'f:l:p:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    fmt = None
    lines = None
    prefix = None
    for o, a in opts:
        if o == '-f':
            if a == 'fasta' or a == 'fastq':
                fmt = a
            else:
                print >> sys.stderr, "Wrong format"
        if o == '-l':
            lines = int(a)
        if o == '-p':
            prefix = a
    if (fmt == None) or (lines == None) or (prefix == None):
        print >> sys.stderr, "Missing options!"
        sys.exit(2)
    line = 0
    seqs = []
    split = 0
    for fin in args:
        for seq in SeqIO.parse(fin, fmt):
            seqs.append(seq)
            line += 1
            if line == lines:
                SeqIO.write(seqs, prefix + (".%d" % split), fmt)
                split += 1
                seqs = []
                line = 0
    if seqs != []:
        SeqIO.write(seqs, prefix + (".%d" % split), fmt)

if __name__ == '__main__':
    main(sys.argv[1:])

