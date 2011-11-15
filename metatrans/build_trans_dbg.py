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

from Bio import SeqIO
import networkx

import trans_dbg

def build_trans_dbg(reads, k):
    """
    Build the de Bruijn graph of transcriptome reads and return the de Bruijn graph
    
    Args:
        reads: FASTQ transcriptome reads
        k: kmer length
    
    Returns:
        The de Bruijn graph formed by the reads using kmers
    """    
    dbg = trans_dbg.TransDBG(networkx.Graph())
    
    for read in reads:  
        for kmer_s in range(len(read) - k + 1):  # kmer_s: kmer start position
            dbg.add_kmer(read, kmer_s, k)
        
    return trans_dbg
