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

'''
count the # of alignments in the output of BWA (SAM)
'''

import HTSeq
import sys

def CountBWASAMAlign(sam_files):
    al_count = 0
    for sam_file in sam_files:
        for align in HTSeq.SAM_Reader(sam_file):
            if align.aligned == True:
                al_count = al_count + 1
    print al_count

if __name__ == '__main__':
    CountBWASAMAlign(sys.argv[1:])
    sys.exit(0)
