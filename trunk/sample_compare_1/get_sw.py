#!/usr/bin/python

import sys
import string

fin = open(sys.argv[1], 'r')

count = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for line in fin.readlines():
    val = float(string.split(line, ' ')[1])
    n = int(string.split(line, ' ')[2])
    if val <= 0.5:
        count[0] = count[0] + n
    elif val <= 1.0:
        count[1] = count[1] + n
    elif val <= 1.5:
        count[2] = count[2] + n
    elif val <= 2.0:
        count[3] = count[3] + n
    elif val <= 2.5:
        count[4] = count[4] + n
    elif val <= 3.0:
        count[5] = count[5] + n
    elif val <= 3.5:
        count[6] = count[6] + n
    elif val <= 4.0:
        count[7] = count[7] + n
    elif val <= 4.5:
        count[8] = count[8] + n
    else: 
        count[9] = count[9] + n

for n in count:
    print n

fin.close()

sys.exit(0)
