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
import shlex
import tempfile
import os
import string
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline

def SWCompare2(f1, f2):
    def StripFasta(fname):
        return ((fname[::-1]).replace("atsaf.", "", 1))[::-1]
        
    results = open(StripFasta(f1) + "." + StripFasta(f2) + ".water", 'w')
    norm_res = open(StripFasta(f1) + "." + StripFasta(f2) + ".normal", 'w')
    f1_tmp = []
    f2_read = list(SeqIO.parse(f2, 'fasta'))
    for seq1 in SeqIO.parse(f1, 'fasta'):
        (tmpf_handle, tmpf_path) = tempfile.mkstemp(prefix = 'sw-2f')  # temp file for storing sequence
        os.close(tmpf_handle)  # close these files
        f1_tmp.append(tmpf_path)
        SeqIO.write([seq1], tmpf_path, 'fasta')
        water_cline = WaterCommandline(asequence = tmpf_path, bsequence = f2
                                       , gapopen = 10, gapextend = 0.5
                                       , outfile = "stdout")
        water_run = subprocess.Popen(shlex.split(water_cline.__str__())
                                     , stdout = subprocess.PIPE)
        f2_i = 0
        for water_out in water_run.stdout.readlines():
            results.write(water_out)
            water_out = string.strip(water_out)
            strip_s = water_out.replace("# Score: ", "", 1)
            if water_out != strip_s:
                # {norm value} {occurrence}
                norm_res.write("%f %d\n" 
                               % (float(strip_s) / min(float(len(seq1.seq)), float(len(f2_read[f2_i].seq)))
                               , int(float(string.split(string.split(f2_read[f2_i].description,' ')[2],'=')[1]) * float(string.split(string.split(seq1.description,' ')[2],'=')[1]))))
                f2_i = f2_i + 1
    results.close()
    norm_res.close()
    for f1_f in f1_tmp:
        os.remove(f1_f)

if __name__ == '__main__':
    SWCompare2(sys.argv[1], sys.argv[2])
    sys.exit(0)

