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
Compares two sets of assembled transcriptome reads (fasta) using SW
'''

'''
TODO: utilize read count
'''

import sys
import subprocess
import tempfile
import os
import shutil
from Bio import SeqIO

def SWCompare2(f1, f2):
    tmpd = tempfile.mkdtemp(prefix = 'sw-2d')  # temp directory for files doing pairwise SW
    f1_tmp = []
    for seq1 in SeqIO.parse(f1, 'fasta'):
        (tmpf_handle, tmpf_path) = tempfile.mkstemp(prefix = 'sw-2f')  # temp file for storing sequence
        f1_tmp.append(tmpf_path)
        SeqIO.write([seq1], tmpf_path, 'fasta')
    shutil.rmtree(tmpd, ignore_errors = True)  # TODO: handle errors
    for f1_f in f1_tmp:
        os.remove(f1_f)

if __name__ == '__main__':
    fout = open("results", 'w')
    SWCompare2(sys.argv[1], sys.argv[2])
    fout.close()
    sys.exit(0)

