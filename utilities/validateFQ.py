#!/usr/bin/env python
#
#
#   validateFQ.py
#   A quickie tool for validating the correctness of a FASTQ file
#
#   Takes an input FASTQ file   
#
#   Usage: python validateFQ.py -i /path/to/FASTQ_file
#
#
# Python 3 ready

import sys
import argparse
from Bio import SeqIO


def func_parser() -> argparse.Namespace:
    """
    Defines what arguments the program requires, and argparse will figure out how to parse those out of sys.argv

    :return: an instance of the argparse class that can be used to access command line arguments
    """

    parser = argparse.ArgumentParser(description='validateFQ.py')
    parser.add_argument('-i', type=str, required=True, metavar='<str>', help="input_file.fastq")
    args = parser.parse_args()

    return args


def main():
    """
    Validates a given FASTQ file for correct format

    :return: None
    """
    
    args = func_parser();
    try:
        for record in SeqIO.parse(args.i, "fastq"):
            pass
        print('\n' + args.i +' Verified!')
    except Exception as e: 
        print(e)
        print('Please use correct format file')


if __name__ == '__main__':
    main()