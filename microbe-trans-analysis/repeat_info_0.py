#!/usr/bin/env python
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

"""
Usage:
./repeat_info_0.py [align sam] 
"""

import sys
import HTSeq
from Bio import SeqIO

def main(argv):
    rep_count = {}
    queries = []
    for align in HTSeq.SAM_Reader(argv[1]):
        if align.aligned:
            if align.read.name not in queries:
                queries.append(align.read.name)
            else:
                if align.read.name not in rep_count:
                    rep_count[align.read.name] = [1, len(align.read.seq)]
                else:
                    rep_count[align.read.name][0] += 1
    print rep_count
    
if __name__ == "__main__":
    main(sys.argv)
    sys.exit(0)

