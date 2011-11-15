#!/usr/bin/python

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
Analysis metatranscriptome sequences using de Bruijn graph
"""

from Bio import SeqIO
import networkx
from matplotlib import pyplot

import read_fastq
import trans_dbg

def dbg_metatrans(k, single_fastq):
    """
    Use de Bruijn graph to analyze metatranscriptome sequences
    
    Args:
        k: kmer length 
        single_fastq: FASTQ files containing single reads
    """   
    single = read_fastq.read_fastq(single_fastq)
    dbg = trans_dbg.TransDBG(single, k)
    
    networkx.draw(dbg.graph)
    pyplot.show()
