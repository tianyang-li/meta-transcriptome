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
read results from blast_relate.py
should be used alone
'''

import sys
import string
from Bio import SeqIO
import json

def GetTransComp(json_file):
    '''
    Return a list of components in dictionary
    '''
    trans_comp = []
    json_dt = open(json_file, 'r')
    for line in json_dt:
        line = string.strip(line)
        comp = json.loads(line)
        trans_comp.append(comp)
    return trans_comp

if __name__ == '__main__':
    print "Cannot be used alone!"
    sys.exit(0)


