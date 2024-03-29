#!/usr/bin/python

#  Copyright (C) 2011 Tianyang Li
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

"""
get ids of seqs that have alignments in db
"""

import sys
from HTSeq import SAM_Reader

def main(args):
    seq_ids=set([])
    for fin in args:
        for align in SAM_Reader(fin):
            if align.aligned:
                seq_ids.add(align.iv.chrom)
    for seq_id in seq_ids:
        print seq_id

if __name__=='__main__':
    main(sys.argv[1:])
