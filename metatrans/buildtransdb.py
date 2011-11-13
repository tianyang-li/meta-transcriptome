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

def BuildTransDB(reads, k):
    """
    Build the de Bruijn graph of transcriptome reads and return the de Bruijn graph
    
    Keyword arguments:
    reads -- FASTQ transcriptome reads
    k -- kmer length
    """

