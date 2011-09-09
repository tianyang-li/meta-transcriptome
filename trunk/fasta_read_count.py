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
fasta_read_count.py {fasta files}
'''

import sys
import string

if __name__ == '__main__':
    read_count = 0
    for fasta in sys.argv[1:]:
        fin = open(fasta, 'r')
        while True:
            line = fin.readline()
            line = string.strip(line)
            if line == '':
                break
            if line[0] == '>':
                read_count = read_count + 1
        fin.close()
    print read_count
    sys.exit(0)


