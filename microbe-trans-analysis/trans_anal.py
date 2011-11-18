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

def trans_anal(trans, seqs):
    # number of transcripts
    print len(seqs)
    trans_comps = networkx.weakly_connected_component_subgraphs(trans.graph)
    # number of connected components in the de Bruijn graph
    print len(trans_comps)
    
    trans_comp_stat = open("trans_comp_stat", 'w')
    
    for trans_comp in trans_comps:
        
        # repeat analysis: 
        #     count the number of strongly connected 
        #     components (which shouldn't occur if no repeats)
        sccs = networkx.strongly_connected_components(trans_comp)
        repeat_count = 0  # not so rigorous???
        for scc in sccs:
            if len(scc) <= 1:
                break
            repeat_count += 1
        
        source_count = 0
        sink_count = 0
        
        for v in trans_comp.nodes():        
            if (list(trans_comp.in_degree_iter([v]))[0][1] == 0):
                source_count += 1
            if (list(trans_comp.out_degree_iter([v]))[0][1] == 0):
                sink_count += 1
        
        # graph size, # of sccs, # of sources, # of sinks
        comp_stat = "%d, %d, %d, %d\n" % (trans_comp.number_of_nodes(), repeat_count, source_count, sink_count)
        trans_comp_stat.write(comp_stat)
        
    trans_comp_stat.close()

def main(argv):
    seqs = read_reads.read_reads(argv[2:], 'fasta')
    trans = trans_dbg.TransDBG(seqs, int(argv[1]))
    trans_anal(trans, seqs)
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)
