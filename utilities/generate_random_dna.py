#!/usr/bin/env python
#
#
#   generate_random_dna.py
#   Generates a random DNA string of given length
#
#   Takes an integer as input and generates a String   
#
#   Usage: python generate_random_dna.py
#
#
# Python 3 ready


import random


def generate_random_dna(lnth: int, seed: int = None) -> str:
    """
    Takes a parameter length and returns a randomly generated DNA string of that length
    :param lnth: how long of a string to generate
    :param seed: Optional seed to produce reproducibly random results
    :return: randomly generated string
    """
    set = ["A", "G", "C", "T"]
    if seed:
        random.seed(seed)
    else:
        random.seed()
    ret = ""
    for i in range(lnth):
        ret += random.choice(set)
    return ret


if __name__ == '__main__':
    print(generate_random_dna(10))
    print(generate_random_dna(10, 1))