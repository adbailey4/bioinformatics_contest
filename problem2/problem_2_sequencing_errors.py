#!/usr/bin/env python
"""Solve sequencing errors """
########################################################################
# File: problem_2_sequencing_errors.py
#  executable: problem_2_sequencing_errors.py
#
# Author: Andrew Bailey
# History: Created 02/04/19
########################################################################

import os
import numpy as np
from scipy.stats import binom

def p_e_given_match_m_c(n_errors, total):
    """Calculate prob_error_given_n_same_errors_and_total_site_coverage p(e | m, c) """
    if total == 0:
        return 1
    else:
        if n_errors == (total / 2.0):
            return 0.5
        elif n_errors > (total / 2.0):
            return 1
        else:
            return 0


def p_m_given_total_prob(m, n, p):
    """Calculate binomial distribution for m errors given n draws with binomial probability of p"""
    rv = binom(n, p)
    return rv.pmf(m)


def p_e_given_n_and_m(m, n, p):
    """Probability of error given number and error prob"""
    return p_e_given_match_m_c(m, n) * p_m_given_total_prob(m, n, p)


def p_e_given_c(n, p):
    """Probability of error given number of aligned bases and probability of error for the base"""
    total_prob = 0
    for m in range(n+1):
        total_prob += p_e_given_n_and_m(m, n, p)

    return total_prob


def p_coverage_given_L_k_n(c, L, n, k):
    """Calculate the probabilty of covering a site in the genome given
    genome length L, kmer length n, and number of reads k"""
    p_coverage = n / L
    rv = binom(k, p_coverage)
    return rv.pmf(c)


def p_error(L, n, p, k):
    """Calculate the probability of error for the genome
    :param L: len genome
    :param n: len reads
    :param k: n reads
    :param p: prob of error
    """
    total_prob = 0
    for c in range(k+1):
        p_c = p_coverage_given_L_k_n(c, L, n, k)
        p_e = p_e_given_c(c, p)
        total_prob += p_c * p_e

    return total_prob * L


def main():
    in_file = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem2/1.txt"
    in_file = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem2/test.txt"

    out_file = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem2/problem2_answers.txt"
    with open(out_file, 'w') as out_f:
        with open(in_file, 'r') as fh:
            iterations = int(fh.readline().rstrip())
            for i in range(iterations):
                data = fh.readline().split()
                error = p_error(int(data[0]), int(data[1]), float(data[2]), int(data[3]))
                print(i, int(data[0]), int(data[1]), float(data[2]), int(data[3]), error)
                # print(error, file=out_f)


if __name__ == '__main__':
    main()
