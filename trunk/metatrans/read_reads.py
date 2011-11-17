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

from Bio import SeqIO
import trans_dbg

def read_reads(in_reads, reads_type): 
    """
    Read reads in a list of FASTQ or FASTA files into memory 
        and return a list containing the reads 
    
    Args:
        reads: list of FASTQ or FASTA files
        reads_type: FASTA or FASTQ
    
    Returns:
        A list containing trans_dbg.Read
    """   
    reads = []
    for reads_file in in_reads:
        for read in SeqIO.parse(reads_file, reads_type):
            reads.append(trans_dbg.Read(read))
    return reads

