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
./repeat_info_0.py [align sam] [usearch uc]
"""

import sys
import HTSeq
import string

def main(argv):
    repeats = {}  
    
    for align in HTSeq.SAM_Reader(argv[1]):
        if align.aligned:
            if align.read.name not in repeats:
                repeats[align.read.name] = [align.iv]
            else:
                repeats[align.read.name].append(align.iv)
                
    uc = open(argv[2], 'r')
    
    clusters = []
    
    for line in uc:
        if line[0] != "#":
            line = string.strip(line)
            fields = line.split("\t")
            if int(fields[1]) not in clusters:
                clusters.append([])
            clusters[int(fields[1])].append(fields[8])
    
    uc.close()
    
    print clusters
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)
