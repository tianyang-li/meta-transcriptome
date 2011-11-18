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
ko_get_genes_fasta.py [file containing KO numbers] [output FASTA file]
"""

import SOAPpy
import sys

def main(argv):
    wsdl = "http://soap.genome.jp/KEGG.wsdl"
    serv = SOAPpy.WSDL.Proxy(wsdl)
    
    ko_list = open(argv[1], 'r')
    ko_genes = open(argv[2], 'w')
    
    count = 0
    
    gene_list = []
    
    for ko_num in ko_list:
        ko_num = ko_num.strip()
        ko_res = serv.get_genes_by_ko("ko:%s" % ko_num, "all")
        for gene_entry in ko_res:
            if gene_entry['entry_id'] not in gene_list:
                gene_list.append(gene_entry['entry_id'])
                fasta_str = serv.bget("-f -n n %s" % gene_entry['entry_id'])
                ko_genes.write(fasta_str)
                
                count += 1
                print "Got gene #%d" % count 
    
    ko_list.close()
    ko_genes.close()
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)



