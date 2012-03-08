#!/usr/bin/env python

# Copyright (C) 2012 Tianyang Li
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

import sys

def main(args):
    for fin in args:
        with open(fin, 'r') as b6:
            for entry in b6:
                entry = entry.strip().split("\t")
                query_id = entry[0].split(" ")[0]
                db_id = entry[1]
                print query_id, db_id
    
if __name__ == '__main__':
    main(sys.argv[1:])









