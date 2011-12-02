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
                repeats[align.read.name] = [align.read.seq, [align.iv]]
            else:
                repeats[align.read.name][1].append(align.iv)
                
    uc = open(argv[2], 'r')
    
    clusters = []
    
    for line in uc:
        if line[0] != "#":
            line = string.strip(line)
            fields = line.split("\t")
            if int(fields[1]) >= len(clusters):
                clusters.append([])
            if fields[8] not in clusters[int(fields[1])]:
                clusters[int(fields[1])].append(fields[8])
    
    uc.close()

    tmp_clusters = []
    for cl in clusters:
        def repeat_in_cluster(my_repeats, my_cl):
            for read in my_cl:
                if read in my_repeats:
                    return read
            return None
        
        clr = repeat_in_cluster(repeats, cl)
        if (len(cl) > 1) and (clr != None):
            reduc_cl = []

            for read in cl:
                if read in repeats:
                    read_good = True
                    for rread in reduc_cl:
                        for iv in repeats[read][1]:
                            for riv in repeats[rread][1]:
                                if (iv.chrom == riv.chrom) and not((iv.start >= riv.end) or (riv.start >= iv.end)):
                                    read_good = False
                    if read_good:
                        reduc_cl.append(read)
            
            if len(reduc_cl) > 1:
                print len(reduc_cl) 
                tmp_clusters.append(reduc_cl)
            
    clusters = tmp_clusters
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)
