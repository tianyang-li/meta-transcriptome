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

from Bio import SeqIO
import sys
import string
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def BlastClassify(fasta_files):
    for fasta_file in fasta_files:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            blast_rec = list(NCBIXML.parse(NCBIWWW.qblast("blastx", "nr", seq.format('fasta'))))
            blast_rec.extend(list(NCBIXML.parse(NCBIWWW.qblast("blastx", "env_nr", seq.format('fasta')))))
            seq_accession = []
            for rec in blast_rec:
                for align in rec.alignments:
                    seq_accession.append(string.split(string.split(align.hit_id, "|")[3], ".")[0])
        
if __name__ == '__main__':
    BlastClassify(sys.argv[1:])
    sys.exit(0)
