#!/usr/bin/env python
#
#
#   compute_gc.py
#   Compute GC and coverage model for gen_reads_runner.py
#
#   Takes output file_list from bedtools genomecov to generate GC/coverage model
#
#   Usage: bedtools genomecov -d -ibam input.bam -g reference.fa > genomeCov.dat
#          python compute_gc.py -r reference.fa -i genomeCov.dat -w [sliding window length] -o output_name
#
#
# Updated to Python 3 standards

import argparse
import gzip
import pickle
import time
import pdb

import numpy as np
from Bio import SeqIO
import pybedtools


def process_fasta(file: str) -> dict:
    """
    Takes a fasta file_list, converts it into a dictionary of upper case sequences. Does some basic error checking,
    like the file_list is readable and the reference dictionary is not empty
    :param file: path to a fasta file_list
    :return: dictionary form of the sequences indexed by chromosome
    """
    ref_dict = {}

    try:
        # reads in fasta file_list, converts sequence to upper case
        ref_dict = {rec.id: rec.seq.upper() for rec in SeqIO.parse(file, "fasta")}
    except UnicodeDecodeError:
        # if the file_list isn't readable, this exception should catch it
        print("Input file_list incorrect: -r should specify the reference fasta")
        exit(1)

    if not ref_dict:
        # if the file_list was readable by SeqIO but wasn't a fasta file_list, this should catch it
        print("Input file_list incorrect: -r should specify the reference fasta")
        exit(1)

    return ref_dict


def process_genomecov(file: str, ref_dict: dict, window: int) -> dict:
    """
    Takes a genomecov file_list and converts it into a dictionary made up of 'window' sized sections
    that record the number of GCs and the coverage measure for each section.
    :param file: path to a genomecov file_list
    :param ref_dict: dictionary created from using the process_fasta function
    :param window: Length of each section of base pairs to count in the reference dictionary
    :return: dictionary form of genomecov file_list based on window size and ref_dict data
    """
    gc_bins = {n: [] for n in range(window + 1)}

    # variables needed to parse coverage file_list
    current_line = 0
    current_ref = None
    current_cov = 0

    f = open(file, 'r')
    for line in f:
        splt = line.strip().split('\t')
        if current_line == 0:
            current_ref = splt[0]
            current_pos = int(splt[1]) - 1

        if current_ref not in ref_dict:
            continue

        current_line += 1
        current_cov += float(splt[2])

        if current_line == window:
            current_line = 0
            seq = str(ref_dict[current_ref][current_pos:current_pos + window])
            if 'N' not in seq:
                gc_count = seq.count('G') + seq.count('C')
                gc_bins[gc_count].append(current_cov)
            current_cov = 0

    f.close()
    return gc_bins


def calculate_coverage(bin_dict: dict, window: int) -> float:
    """
    Takes the dictionary created in process_genomecov and finds the average coverage value.
    Also ouputs the average coverage value for each window, along with the number of entries in that window.
    :param bin_dict: dictionary created from using the process_genomecov function
    :param window: Length of each section of base pairs to count, 
                   should be the same as the window value in process_genomecov
    :return: Average coverage value for the whole sample, along with average coverage values for each window.
    """
    running_total = 0
    all_mean = 0.0
    for k in sorted(bin_dict.keys()):
        if len(bin_dict[k]) == 0:
            print('{0:0.2%}'.format(k / float(window)), 0.0, 0)
            bin_dict[k] = 0
        else:
            my_mean = np.mean(bin_dict[k])
            my_len = len(bin_dict[k])
            print('{0:0.2%}'.format(k / float(window)), my_mean, my_len)
            all_mean += my_mean * my_len
            running_total += my_len
            bin_dict[k] = my_mean

    return all_mean / float(running_total)


def main():
    """
    Reads in arguments and processes the inputs to a GC count for the sequence.
    Parameters:
        -i is the genome coverage input file_list (genomecov)
        -r is the reference file_list (fasta)
        -o is the prefix for the output
        -w is the sliding window length. The default is 50, but you can declare any reasonable integer
    :return: None
    """
    parser = argparse.ArgumentParser(description='compute_gc.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i', type=str, required=True, metavar='input', help="input.genomecov")
    parser.add_argument('-r', type=str, required=True, metavar='reference', help="reference.fasta")
    parser.add_argument('-o', type=str, required=True, metavar='output prefix',
                        help="prefix for output (/path/to/output)")
    parser.add_argument('-w', type=int, required=False, metavar='sliding window',
                        help="sliding window length [50]", default=50)
    args = parser.parse_args()

    (in_gcb, ref_file, window_size, out_p) = (args.i, args.r, args.w, args.o)

    print('Reading ref...')
    allrefs = process_fasta(ref_file)

    tt = time.time()
    print('Reading genome coverage file_list...')
    gc_bins = process_genomecov(in_gcb, allrefs, window_size)

    print("Calculating average coverage...")
    average_coverage = calculate_coverage(gc_bins, window_size)

    print('AVERAGE COVERAGE =', average_coverage)

    y_out = []
    for k in sorted(gc_bins.keys()):
        gc_bins[k] /= average_coverage
        y_out.append(gc_bins[k])

    print('saving model...')
    pickle.dump([range(window_size + 1), y_out], gzip.open(out_p + ".pickle.gz", 'wb'))
    print(time.time() - tt, '(sec)')


if __name__ == "__main__":
    main()
