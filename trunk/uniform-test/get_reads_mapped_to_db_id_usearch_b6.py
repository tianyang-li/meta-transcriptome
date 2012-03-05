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

def main(args):
    fmt, db, b6 = None, None, None, None
    try:
        opts, reads_fins = getopt.getopt(args, 'f:d:b:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-f':
            # format of reads
            fmt = a
        if o == '-d':
            # file containing all the db seqs, assumed to be in FASTA
            db = a
        if o == '-b':
            # usearch b6 files saying which reads mapped to db seqs
            b6 = a
    if fmt == None or db == None or b6 == None or reads_fins == []:
        print >> sys.stderr, "missing arguments"
        sys.exit(2)



if __name__ == '__main__':
    main(sys.argv[1:])

