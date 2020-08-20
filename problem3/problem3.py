#!/usr/bin/env python
"""Solve problem3 """
########################################################################
# File: problem3.py
#  executable: problem3.py
#
# Author: Andrew Bailey
# History: Created 02/08/19
########################################################################

from py3helpers.utils import find_substring_indices
from py3helpers.seq_tools import get_pairwise_cigar, pairwise_alignment_accuracy
import operator
stats = {'a':1000, 'b':3000, 'c': 100, 'd':3000}


def get_kmers_from_seq(sequence, k):
    list_of_kmers = []
    for x in range(len(sequence)-k+1):
        list_of_kmers.append(sequence[x:k+x])
    return list_of_kmers

def get_nonoverlapping_kmers(sequence, k):
    list_of_kmers = []
    for i in range(k):
        non_overlap = []
        for x in range(i, len(sequence)-k+1, k):
            non_overlap.append(sequence[x:k+x])
        if non_overlap != []:
            list_of_kmers.append(non_overlap)
    return list_of_kmers


def find_perfect_match_kmers(list_of_kmers):
    kmer_dict = {kmer: 0 for kmer in list_of_kmers}
    for kmer in set(list_of_kmers):
        for kmer2 in list_of_kmers:
            if kmer == kmer2:
                kmer_dict[kmer] += 1
    return kmer_dict


def find_answer(test_seq, n, l, d):
    kmers = get_kmers_from_seq(test_seq, k=l)
    kmer_dict = {kmer: 0 for kmer in kmers}
    for kmer in set(kmers):
        for kmer2 in kmers:
            pairwise_alignment_accuracy(kmer, kmer2)
                kmer_dict[kmer] += 1
    return kmer_dict

    return


def write_file(test_seq, kmer, outfile):
    with open(outfile, "w") as fh:
        indices = [x+1 for x in find_substring_indices(test_seq, kmer, overlap=False)]
        fh.write(kmer+"\n")
        for x in indices:
            fh.write(str(x)+" ")
            fh.write(str(len(kmer)))
            fh.write("M"+"\n")

def main():
    test1 = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem3/2.txt"
    with open(test1, 'r') as fh:
        params = [int(x) for x in fh.readline().split()]
        n = params[0]
        l = params[1]
        d = params[2]
        print(n, l, d)
        sequence = fh.readline()

    test_seq = "AATGGGACACATGCGCTGGGAGCCTGGTAATAAGCTGATTGAACTACAGATGACCCGCAAATGGAGACCTTTAGGAAAGAGTATCAAGGAAGTTAGGCGACACACGTACGAAGTGCGCCCAGATCTGACTTAAGAAACGTCGGGGTCATTTGGATACTAAGTCAAGCGAGAGCACGACACCCGCATTCGACCAGTGACCGAATGGGACACATGCGCTGGGAGCCTGCGACGTTCGCCGGCGGTAACGGCTTAACGGGGCTTGTTGCTGCTAGTCGGCGATATAGGTCTTCAGTAAAGCCATCTACTGGCCGCTTTTGAATGGTACCGAAGAGCAAAGCAAGTTCATTTGATTATTCTACTGTGGCGATTTCTATTCGTCGTGTTATAACATTGATTGCTCGCGATCGGGCCCGTTAGGCTTACTTCTGCGGAACGTGTTCTCGAAGGGATAGGTGCGAGGTGCGGGGCATGGAATTTTAGTCCTCCCTCTCCAAGCTGGCCGCTCTTCATGTTGTCATTTTTAGAATTTGGGTTGAGGTCCCCGCATAACAAACACTTTGGGACACATGCGCTGGGAGCCTCCAAAGGCAGTAGGCTTGGGGCACATGGGACACATGCGCTGGGAGCCTACAGGAACCGCTCTCACACGGTCCCGAAATTTGCCCGTGTGACCAACAACATCCTTTTATTTGTCGGCTGAAGTCATTGGTAGCGTGTTCACCCTTACGTGCGTAACCCCAGCGCGAATCTTCACCCCTAATAGTGCCGAGTACAGCTGGGTACCCGCTGCCAGATGAATGTACTAAGTCGGAAGGCATCTGTTTATCTGAGGAGCATTGCCTGCGGGCATAAAAATGGGACACATGCGCTGGGAGCCTCTGTTAGTTGCGAACTACGGACATGGTCCGACACTAGAAAGATTTGTATGGAAGCGGATCGAAGCCCCTGCTTCGACTGTACACCCCATGTCCCGTTCTGAACGATGGGACACATGCGCTGGGAGCCTCCGGGCAATTATGACCACAACTTCGGAGGGTTGTAGATCGGATTATTGCGTATCGTCGCAGTTTTTCCACACGGACTATCGTCGTCTAAACTAACCGGGGGGGCTCAGGTGGGACACATGCGCTGGGAGCCTCTG"
    outfile = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem3/2answer.txt"
    write_file(test_seq, "TGGGACACATGCGCTGGGAGCCT", outfile)

if __name__ == '__main__':
    main()
