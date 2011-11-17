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

import sys
import networkx

import read_reads
import trans_dbg

def main(argv):
    seqs = read_reads.read_reads(argv[2:], 'fasta')
    trans = trans_dbg.TransDBG(seqs, int(argv[1]))
    
    # number of transcripts
    print len(seqs),
    trans_comp = networkx.weakly_connected_component_subgraphs(trans.graph)
    # number of connected components in the de Bruijn graph
    print len(trans_comp)
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)
