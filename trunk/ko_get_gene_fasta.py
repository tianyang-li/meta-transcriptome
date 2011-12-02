#!/data1/tyli/bin/python2.6

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

from suds import client
import sys
import os

def main(argv):
    wsdl = "http://soap.genome.jp/KEGG.wsdl"
    proxyOpts = dict(http="bio_user:12345678@166.111.74.4:3128", https="bio_user:12345678@166.111.74.4:3128") 
    serv = client.Client(wsdl, proxy=proxyOpts)

    ko_list = open(argv[1], 'r')
    ko_genes = open(argv[2], 'w')
    
    count = 0
    
    gene_list = []
    
    for ko_num in ko_list:
        ko_num = ko_num.strip()
        while True:
            try:
                print "Getting KO gene list: %s" % ko_num
                ko_res = serv.service.get_genes_by_ko("ko:%s" % ko_num, "all")
                break
            except BaseException as err:
                print >> sys.stderr, str(err)
        for gene_entry in ko_res:
            if str(gene_entry['entry_id']) not in gene_list:
                gene_list.append(str(gene_entry['entry_id']))
                while True:
                    try:
                        print "Getting gene: %s in %s" % (str(gene_entry['entry_id']), ko_num)
                        fasta_str = serv.service.bget("-f -n n %s" % str(gene_entry['entry_id']))
                        break
                    except BaseException as err:
                        print >> sys.stderr, str(err)
                ko_genes.write(str(fasta_str))
                
                count += 1
                print "Got gene #%d" % count 
    
    ko_list.close()
    ko_genes.close()
    
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)



