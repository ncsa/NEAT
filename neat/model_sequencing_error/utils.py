"""
Utilities to generate the sequencing error model
"""

import logging
import numpy as np
import argparse
import sys
import pickle
import pathlib
import pysam
from functools import reduce


from bisect import bisect_left

__all__ = [
    "take_closest",
    "parse_file"
]

_LOG = logging.getLogger(__name__)


def take_closest(bins, quality):
    """
    Assumes bins is sorted. Returns closest value to quality.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(bins, quality)
    if pos == 0:
        return bins[0]
    if pos == len(bins):
        return bins[-1]
    before = bins[pos - 1]
    after = bins[pos]
    if after - quality < quality - before:
        return after
    else:
        return before


def parse_file(input_file, quality_scores, off_q, max_reads, n_samp):
    _LOG.info("Hello World!")
    # # Takes a gzip or sam file and returns the simulation's average error rate,
    # print('reading ' + input_file + '...')
    # is_aligned = False
    # lines_to_read = 0
    # try:
    #     if input_file[-4:] == '.bam' or input_file[-4:] == '.sam':
    #         print('detected aligned file....')
    #         stats = pysam.idxstats(input_file).strip().split('\n')
    #         lines_to_read = reduce(lambda x, y: x + y, [eval('+'.join(l.rstrip('\n').split('\t')[2:])) for l in stats])
    #         f = pysam.AlignmentFile(input_file)
    #         is_aligned = True
    #     else:
    #         print('detected fastq file....')
    #         with pysam.FastxFile(input_file) as f:
    #             for _ in f:
    #                 lines_to_read += 1
    #         f = pysam.FastxFile(input_file)
    # except FileNotFoundError:
    #     print("Check input file. Must be fastq, gzipped fastq, or bam/sam file.")
    #     sys.exit(1)
    #
    # actual_readlen = 0
    # current_line = 0
    # quarters = lines_to_read // 4
    #
    # # A list of n dictionaries, n being the length of a sequence // Put outside this loop
    # error_model = {
    #     # a list of q scores
    #     'quality_scores': quality_scores,
    #     # number of bins * length of sequence data frame to hold the quality score probabilities
    #     'quality_score_probabilities': [],
    #     'quality_offset': off_q,
    #     'avg_error': {},
    #     'error_parameters': [
    #         # sequencing substitution transition probabilities
    #         [[0., 0.4918, 0.3377, 0.1705], [0.5238, 0., 0.2661, 0.2101], [0.3754, 0.2355, 0., 0.3890],
    #          [0.2505, 0.2552, 0.4942, 0.]],
    #         # if a sequencing error occurs, what are the odds it's an indel?
    #         0.01,
    #         # sequencing indel error length distribution
    #         [0.999, 0.001],
    #         [1, 2],
    #         # Given an indel error occurs, what are the odds it's an insertion?
    #         0.4,
    #         # Given an insertion error occurs, what's the probability of it being an A, C, G, T?
    #         [0.25, 0.25, 0.25, 0.25]
    #     ]
    # }
    #
    # if is_aligned:
    #     g = f.fetch()
    # else:
    #     g = f
    #
    # # Used to get read length, the read length will = the fist read length.
    # obtained_read_length = False
    # temp_q_count = 0
    #
    # for read in g:
    #     if is_aligned:
    #         qualities_to_check = read.query_alignment_qualities
    #     else:
    #         qualities_to_check = read.get_quality_array()
    #
    #     # read length = the length of the first read
    #     if actual_readlen == 0:
    #         actual_readlen = len(qualities_to_check) - 1
    #         print('assuming read length is uniform...')
    #         print('detected read length (from first read found):', actual_readlen)
    #
    #     # check if read length is more than 0 and if we have the read length already**
    #     if actual_readlen > 0 and not obtained_read_length:
    #         temp_q_count = np.zeros((actual_readlen, len(quality_scores)))
    #         error_model['quality_score_probabilities'] = np.zeros((actual_readlen, len(quality_scores)), dtype=float)
    #         obtained_read_length = True
    #
    #     # sanity-check readlengths
    #     if len(qualities_to_check) - 1 != actual_readlen:
    #         print('skipping read with unexpected length...')
    #         continue
    #
    #     for i in range(0, actual_readlen):
    #         # The qualities of each base
    #         q = qualities_to_check[i]
    #         bin = take_closest(quality_scores, q)
    #         bin_index = quality_scores.index(bin)
    #         temp_q_count[i][bin_index] += 1
    #
    #     # loading
    #     current_line += 1
    #     if current_line % quarters == 0:
    #         print(f'{(current_line / lines_to_read) * 100:.0f}%')
    #     if 0 < max_reads <= current_line:
    #         break
    # f.close()
    #
    # # Probability calculator
    # base_index = 0
    # # for every dictionary(scores of a single base in the reads)
    # for base in temp_q_count:
    #     total = sum(base)
    #     bin_index = 0
    #     for score in base:
    #         error_model['quality_score_probabilities'][base_index][bin_index] = score / total
    #         bin_index += 1
    #     base_index += 1
    #
    # # A Discrete distribution of the quality_score probabilities
    # Discretes = []
    # for base in error_model['quality_score_probabilities']:
    #     Discretes.append(DiscreteDistribution(error_model['quality_scores'], base))
    #
    # # A counter for the Discrete distribution results run for n_samp times
    # count_dict = {}
    # for q in error_model['quality_scores']:
    #     count_dict[q] = 0
    # # example: {0: 0, 12: 0, 24: 0, 36: 0}
    # lines_to_sample = len(range(1, n_samp + 1))
    # # Divide the reads into 1/4s for the loading bar
    # samp_quarters = lines_to_sample // 4
    # for samp in range(1, n_samp + 1):
    #     if samp % samp_quarters == 0:
    #         # loading bar
    #         print(f'{(samp / lines_to_sample) * 100:.0f}%')
    #     for i in range(actual_readlen):
    #         my_q = Discretes[i].sample()
    #         my_q = take_closest(quality_scores, my_q)
    #         count_dict[my_q] += 1
    #
    # print(count_dict)
    #
    # # Calculates the average error rate
    # tot_bases = float(sum(count_dict.values()))
    # avg_err = 0.
    # for k in sorted(count_dict.keys()):
    #     eVal = 10. ** (-k / 10.)
    #     print(k, eVal, count_dict[k])
    #     avg_err += eVal * (count_dict[k] / tot_bases)
    # print('AVG ERROR RATE:', avg_err)
    #
    # return error_model