#!/usr/bin/env python

#  Copyright (C) 2012 Tianyang Li
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
get db seqs that are completely covered by query seqs
"""

import getopt
import sys

def main(args):
    query_ids, db_ids, pair_ids = None, None, None
    try:
        opts, b6s = getopt.getopt(args, 'q:d:p:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-q':
            # write query ids to this file
            query_ids = arg
        if opt == '-d':
            # write db ids to this file
            db_ids = arg
        if opt == '-p':
            # write 
            # query_id db_id
            # to this file
            pair_ids = arg
    if pair_ids == None or query_ids == None or db_ids == None or b6s == []:
        print >> sys.stderr, "Missing arguments"
        sys.exit(2)
        
    query_ids = open(query_ids, 'r')
    db_ids = open(db_ids, 'r')
    pair_ids = open(pair_ids, 'r')
    
    for fin in b6s:
        with open(fin, 'r') as b6:
            for entry in b6:
                entry = entry.strip().split("\t")
                #TODO
                
    query_ids.close()
    db_ids.close()
    pair_ids.close()
    
    
if __name__ == '__main__':
    main(sys.argv[1:])








