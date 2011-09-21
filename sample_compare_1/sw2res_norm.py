#!/usr/bin/python

'''
SW scores from the output of sw_2.py is normalized in the following way:
    SW_score / align_len
'''

import sys
import string

def SWNorm(swres):
    def StripWater(fname):
        return ((fname[::-1]).replace("retaw.", "", 1))[::-1]
        
    fout = open(StripWater(swres) + ".norm", 'w')
    fin = open(swres, 'r')
    
    while True:
        line = fin.readline()
        if line == "":
            break
        line = string.strip(line)
        strip_len = line.replace("# Length: ", "", 1)
        if line != strip_len:
            align_len = float(strip_len)
            fin.readline()  # identity
            fin.readline()  # similarity
            fin.readline()  # gaps
            line = fin.readline()
            line = string.strip(line)
            line = line.replace("# Score: ", "", 1)
            sw_score = float(line)
            norm_s = sw_score / align_len
            fout.write("%f\n" % norm_s)
        
    fin.close()
    fout.close()

if __name__ == '__main__':
    SWNorm(sys.argv[1])
    sys.exit(0)