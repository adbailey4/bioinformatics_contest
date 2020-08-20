#!/usr/bin/env python
"""Solve bee problem"""
########################################################################
# File: problem_1_bee.py
#  executable: problem_1_bee.py
#
# Author: Andrew Bailey
# History: Created 02/04/19
########################################################################

import os
import numpy as np


def next_day_num_bees(n_i, a, b):
    """Calculate the number of bees on the next day"""
    return (n_i*a) - (b*(n_i**2))

def get_limit(n_i, a, b):
    """Get the limit"""
    counter = 0
    running_average = []
    numbers = []
    previous = None
    difference = 100
    while difference > 1e-15:
        counter += 1
        previous = n_i
        n_i = next_day_num_bees(n_i, a, b)
        if n_i < 0:
            return 0
        if n_i > 10e10:
            return -1
        numbers.append(n_i)
        difference = np.abs(n_i-previous)
        running_average.append(difference)
        if counter == 100:
            # print(numbers)
            if np.std(running_average) < 0.000000001:
                if min(numbers) == numbers[-1] and max(numbers) == numbers[0]:
                    return 0
                return np.round(np.mean(numbers[:10]), decimals=5)
            running_average = []
            numbers = []
            counter = 0
    return np.round(n_i, decimals=5)


def main():
    in_file = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem1/tests.txt"
    # in_file = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem1/test1.txt"

    out_file = "/Users/andrewbailey/CLionProjects/bioinformatics_contest/problem1/problem1_answers.txt"
    with open(out_file, 'w') as out_f:
        with open(in_file, 'r') as fh:
            iterations = int(fh.readline().rstrip())
            print(iterations)
            for i in range(iterations):
                data = fh.readline().split()
                limit = get_limit(float(data[0]), float(data[1]), float(data[2]))
                print(i, float(data[0]), float(data[1]), float(data[2]), limit)
                print(limit, file=out_f)

if __name__ == '__main__':
    main()