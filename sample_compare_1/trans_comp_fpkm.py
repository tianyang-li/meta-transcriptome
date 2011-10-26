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

from Bio import SeqIO
import sys

import blast_comp_read

def AnalFPKM(json_file):
    trans_comp = blast_comp_read.GetTransComp(json_file)
    for comp in trans_comp:
        print 1

if __name__ == '__main__':
    AnalFPKM(sys.argv[1])
    sys.exit(0)

