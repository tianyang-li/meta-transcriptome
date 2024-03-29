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
sw_2.py {fasta file 1} {fasta file 2} {similar cutoff}
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
import random
import numpy
import math
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline

def SWCompare2(f1, f2, cutoff):
    K = 0.620391 
    Lambda = 1.33249
    
    def StripFasta(fname):
        return ((fname[::-1]).replace("atsaf.", "", 1))[::-1]
        
    file_name = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(20))  # file names for a run
    results = open(file_name + ".water", 'w')
    norm_res = open(file_name + ".normal", 'w')
    f1_tmp = []
    f2_read = list(SeqIO.parse(f2, 'fasta'))
    d_pval = file_name + "-pval"
    os.mkdir(d_pval)
    f1_i = 0
    for seq1 in SeqIO.parse(f1, 'fasta'):
        (tmpf_handle, tmpf_path) = tempfile.mkstemp(prefix='sw-2f')  # temp file for storing sequence
        os.close(tmpf_handle)  # close these files
        f1_tmp.append(tmpf_path)
        SeqIO.write([seq1], tmpf_path, 'fasta')
        water_cline = WaterCommandline(asequence=tmpf_path, bsequence=f2
                                       , gapopen=5, gapextend=2
                                       , outfile="stdout")
        water_run = subprocess.Popen(shlex.split(water_cline.__str__() + r" -datafile ../EMBOSS-6.4.0/emboss/data/EDNAFULL")
                                     , stdout=subprocess.PIPE)
        f2_i = 0
        al_len = 0
        for water_out in water_run.stdout.readlines():
            results.write(water_out)
            water_out = string.strip(water_out)
            strip_len = water_out.replace("# Length: ", "", 1)
            if water_out != strip_len:
                al_len = int(strip_len)
            strip_s = water_out.replace("# Score: ", "", 1)
            if water_out != strip_s:
                al_score = float(strip_s)
                # {norm - SW/min(a,b)} {norm - SW/al_len} {occurrence}
                norm1 = al_score / min(float(len(seq1.seq)), float(len(f2_read[f2_i].seq)))
                norm2 = al_score / float(al_len)
                pval = 0.0 - numpy.expm1(0.0 - K * len(seq1.seq) * len(f2_read[f2_i].seq) * math.exp(0.0 - Lambda * al_score))
                norm_res.write("%f %f %f %d\n" 
                               % (norm1, norm2, pval
                               , int(float(string.split(string.split(f2_read[f2_i].description, " ")[2], "=")[1]) * float(string.split(string.split(seq1.description, " ")[2], "=")[1]))))
                if pval <= cutoff:
                    f_sim_seq = d_pval + "/" + ("%d-%d" % (f1_i, f2_i))
                    SeqIO.write([seq1], f_sim_seq + "-a", 'fasta')
                    SeqIO.write([f2_read[f2_i]], f_sim_seq + "-b", 'fasta')
                f2_i = f2_i + 1
        f1_i = f1_i + 1
    results.close()
    norm_res.close()
    for f1_f in f1_tmp:
        os.remove(f1_f)

if __name__ == '__main__':
    SWCompare2(sys.argv[1], sys.argv[2], float(sys.argv[3]))
    sys.exit(0)

