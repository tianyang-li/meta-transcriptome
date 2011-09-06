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
fastq2fasta.py {*.fastq files} 
'''

import sys
import fileinput
import string

if __name__ == '__main__':
    for fastq_name in sys.argv[1:]:
        fastq = open(fastq_name, 'r')
        fasta = open(((fastq_name[::-1]).replace('qtsaf', 'atsaf', 1))[::-1], 'w')
        
        state = 0
        while True:
            line = fastq.readline()
            line = string.strip(line)
            if line == '':
                break
            if state == 0:
                line = '>' + line[1:]
            if state == 1 or state == 0:
                fasta.write("%s\n" % (line))
            state = (state + 1) % 4

        fasta.close()
        fastq.close()

    sys.exit(0)

