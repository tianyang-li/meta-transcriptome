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
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
import networkx

class Read(object):
    """
    FASTQ or FASTA read
    
    Attributes:
        fastq: a FASTQ or FASTA read (SeqRecord)
    """
    
    def __init__(self, read):
        self.read = read

class TransDBGNode(object):
    """
    Node in the transcriptome de Bruijn graph that contains a kmer
    
    Attributes: 
        reads: a list of trans_dbg.Read that kmer came from, 
                each entry in reads is a list of the following format:
                    [read, kmer_s]
                read is the read where the kmer came from, 
                kmer_s is the kmer's start position
    """
    
    def __init__(self):
        self.reads = []
        
class TransDBG(object):
    """
    Transcriptome de Bruijn graph
    
    Attributes:
        graph: networkx.DiGraph
                edge direction in the graph is determined as follows:                       
                    ACTGACTGACTG      kmer_node_1
                       ------->
                     CTGACTGACTGA     kmer_node_2
        kmer_dict: maps a kmer to its node in graph
    """
    
    def __init__(self, reads, k): 
        self.graph = networkx.DiGraph()
        self.kmer_dict = {}
        self._build_dbg(reads, k)
    
    def _build_dbg(self, reads, k):
        self._build_dbg_nodes(reads, k)
        self._build_dbg_edges()
    
    def _build_dbg_nodes(self, reads, k):
        for read in reads:
            for kmer_s in range(len(read.read) - k + 1):  # kmer_s: kmer start position
                kmer = str(read.read.seq)[kmer_s : kmer_s + k - 1]
                if kmer not in self.kmer_dict:
                    kmer_node = TransDBGNode()
                    kmer_node.reads.append([read, kmer_s])
                    self.kmer_dict[kmer] = kmer_node
                    self.graph.add_node(kmer_node, kmer=kmer)
                
    def _build_dbg_edges(self):
        for kmer in self.kmer_dict:
            for nuc in IUPACUnambiguousDNA.letters:
                kmer_next = kmer[1:] + nuc
                if kmer_next in self.kmer_dict:
                    self.graph.add_edge(self.kmer_dict[kmer], self.kmer_dict[kmer_next])
