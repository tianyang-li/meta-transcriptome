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
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Author Tianyang Li
#  E-Mail tmy1018 gmail com

'''
Compares two sets of assembled transcriptome reads using SW
'''

'''
TODO: utilize read count
'''

import sys

def SWCompare2(reads1, reads2):
    for read1 in reads1:
        for read2 in reads2:
            # compare read1, read2

if __name__ == '__main__':
    f1 = open(sys.argv[1], 'r')
    f2 = open(sys.argv[2], 'r')
    fout = open("results", 'w')

    f1.close()
    f2.close()
    fout.close()
    sys.exit(0)

