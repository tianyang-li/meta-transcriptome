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
Analyze kmer frequencies and repeats in FASTA or FASTQ sequences

Usage:
./seq_kmer.py [read type: fasta/fastq] [kmer length] 
    [read file 1 (fasta or fastq)] [read file 2 (fasta or fastq)] ...
"""

from Bio import SeqIO
import sys

def SingleSeqKmer(read_str, k):
    """
    Kmer statistics for read_str
    
    Arguments:
        read_str: a string of ACGT sequences
        
    Returns:
        dictionary containing kmer frequencies
    """
    kmer_count = {}
    
    for kmer_s in range(len(read_str) - k + 1):
        kmer = read_str[kmer_s : kmer_s + k - 1]
        if kmer not in kmer_count:
            kmer_count[kmer] = 1
        else:
            kmer_count[kmer] += 1
        
    return kmer_count

def GetKmers(read_str, k):
    """
    Get a list of kmers in read_str
    
    Arguments:
        read_str: a string of ACGT sequences
        
    Returns:
        a list containing all the kmers in read_str
    """
    kmers = []
    for kmer_s in range(len(read_str) - k + 1):
        kmer = read_str[kmer_s : kmer_s + k - 1]
        kmers.append(kmer)
    return kmers

def main(argv):
    k = int(argv[2])
    
    # a list of [read, 
    #     , kmer_self_repeat (int), read_overlap_kmer (int)]  
    reads = []  
    
    for rf in argv[3:]:
        for read in SeqIO.parse(rf, argv[1]):
            reads.append([read])
    
    for read_entry in reads:
        kmer_count = SingleSeqKmer(str(read_entry[0].seq), k)
        
        kmer_self_repeat = 0
        for kmer, occur in kmer_count.iteritems():
            if occur > 1:
                kmer_self_repeat += 1
        read_entry.append(kmer_self_repeat)
        
        read_overlap_kmer = 0
        read_entry.append(read_overlap_kmer)
    
    for ent1 in range(len(reads)):
        for ent2 in range(ent1 + 1, len(reads)):
            set1 = set(GetKmers(str(reads[ent1][0].seq), k))
            set2 = set(GetKmers(str(reads[ent2][0].seq), k))
            if len(set1.intersection(set2)) != 0:
                reads[ent1][2] += 1
                reads[ent2][2] += 1  
    
    for read_entry in reads:
        print read_entry[1], read_entry[2]
            
        
if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)


