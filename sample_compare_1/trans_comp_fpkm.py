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

import sys
import string
import numpy

import blast_comp_read

def AnalFPKM(json_file):
    '''
    return a list of lists of fpkm in a single component
    '''
    trans_comp = blast_comp_read.GetTransComp(json_file)
    trans_comp_fpkm = []
    for comp in trans_comp:
        comp_fpkm = []
        for  i in range(len(comp['trans'])):
            fpkm = float(string.split(string.split(string.split(comp['trans']['%d' % i]['seq'], "\n")[0], " ")[2], "=")[1])
            comp_fpkm.append(fpkm)
        trans_comp_fpkm.append(comp_fpkm)
        if len(comp_fpkm) > 3:
            '''
            mean std_dev std_dev/mean
            '''    
            print "%f %f %f" % (numpy.average(comp_fpkm)
                                , numpy.std(comp_fpkm)
                                , numpy.std(comp_fpkm) / numpy.average(comp_fpkm))
    return trans_comp_fpkm

if __name__ == '__main__':
    AnalFPKM(sys.argv[1])
    sys.exit(0)

