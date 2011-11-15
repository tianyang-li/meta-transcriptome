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

class Read(object):
    """
    FASTQ read
    
    Attributes:
        fastq: a FASTQ read
    """
    
    def __init__(self, fastq):
        self.fastq = fastq

class TransDBGNode(object):
    """
    Node in the transcriptome de Bruijn graph that contains a kmer
    
    Attributes: 
        kmer:
        reads: a list of trans_dbg.Read that kmer came from, 
                each entry in reads is a list of the following format:
                    [read, kmer_s]
                read is the read where the kmer came from, 
                kmer_s is the kmer's start position
    """
    
    def __init__(self, kmer):
        self.kmer = kmer
        self.reads = []
        
class TransDBG(object):
    """
    Transcriptome de Bruijn graph
    
    Attributes:
        graph: networkx.Graph
        kmer_dict: maps a kmer to its node in graph
    """
    
    def __init__(self, G): 
        self.graph = G  
        self.kmer_dict = {}
    
    def add_kmer(self, read, kmer_s, k):
        kmer = str(read.seq)[kmer_s : kmer_s + k - 1]
        kmer_node = TransDBGNode(kmer)
        kmer_node.reads.append([read, kmer_s])
        self.kmer_dict[kmer] = kmer_node
        self.graph.add_node(kmer_node)
