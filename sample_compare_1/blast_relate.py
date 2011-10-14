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

'''
relate each file containing transcript sequences to each other using blastx
'''

from Bio import SeqIO
import sys
import string
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from networkx import nx
from networkx.algorithms.components.connected import number_connected_components

def BlastClassify(fasta_files):
    ac = []  # a list of accessions
    ac_gr = nx.Graph()  # graph where each node represents accessions within a single file 
    for fasta_file in fasta_files:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            blast_rec = list(NCBIXML.parse(NCBIWWW.qblast("blastx", "nr", seq.format('fasta'))))
            blast_rec.extend(list(NCBIXML.parse(NCBIWWW.qblast("blastx", "env_nr", seq.format('fasta')))))
            seq_accession = []
            for rec in blast_rec:
                for align in rec.alignments:
                    seq_accession.append(string.split(string.split(align.hit_id, "|")[3], ".")[0])
            if seq_accession != []:
                ac_gr.add_node(len(ac))
                ac.append(seq_accession)
    for ac1 in ac_gr.nodes():
        for ac2 in ac_gr.nodes():
            if ac1 != ac2:
                if len(set(ac(ac1)) & set(ac(ac2))) != 0:
                    ac_gr.add_edge(ac1, ac2)
    print number_connected_components(ac_gr)
        
if __name__ == '__main__':
    BlastClassify(sys.argv[1:])
    sys.exit(0)
