"""
Utilities to generate the sequencing error model
"""

import logging
import numpy as np

from Bio import SeqIO
from bisect import bisect_left
from scipy.stats import mode
from ..common import open_input

__all__ = [
    "parse_file"
]

_LOG = logging.getLogger(__name__)


def take_closest(bins, quality):
    """
    Assumes bins is sorted. Returns the closest value to quality.

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


def convert_quality_string(qual_str: str, offset: int):
    """
    Converts a plain quality string to a list of numerical equivalents

    :param qual_str: The string to convert
    :param offset: the quality offset for conversion for this fastq
    :return list: a list of numeric quality scores
    """
    ret_list = []
    for i in range(len(qual_str)):
        ret_list.append(ord(qual_str[i]) - offset)

    return ret_list


def expand_counts(count_array: list, scores: list):
    """
    Expands a counting list out into the full tally

    :param count_array: the list to expand
    :param scores: The factors by which to expand the list
    :return np.ndarray: a one-dimensional array reflecting the expanded count
    """
    if len(count_array) != len(scores):
        raise ValueError("Count array and scores have different lengths.")

    ret_list = []
    for i in range(len(count_array)):
        ret_list.extend([scores[i]] * count_array[i])

    return np.array(ret_list)


def parse_file(input_file: str, quality_scores: list, max_reads: int, qual_offset: int):
    """
    Parses an individual file for statistics

    :param input_file: The input file to process
    :param quality_scores: A list of potential quality scores
    :param max_reads: Max number of reads to process for this file
    :param qual_offset: The offset score for this fastq file. We assume the Illumina default of 33.
    :return:
    """

    _LOG.info(f'reading {input_file}')

    fastq_index = SeqIO.index(input_file, 'fastq')

    readlens = []

    counter = 0
    for read_name in fastq_index:
        read = fastq_index[read_name]
        if read.letter_annotations['phred_quality']:
            readlens.append(len(read))
            counter += 1
            if counter >= 1000:
                # takes too long and uses too much memory to read all of them, so let's just get a sample.
                break

    readlens = np.array(readlens)

    # Using the statistical mode seems like the right approach here. We expect the readlens to be roughly the same.
    readlen_mode = mode(readlens, axis=None, keepdims=False)
    if int(readlen_mode.count) < (0.5 * len(readlens)):
        _LOG.warning("Highly variable read lengths detected. Results may be less than ideal.")
    if int(readlen_mode.count) < 20:
        raise ValueError(f"Dataset is too scarce or inconsistent to make a model. Try a different input.")

    read_length = int(readlen_mode.mode)

    _LOG.debug(f'Read len of {read_length}, over {counter} samples')

    total_records_to_read = min(len(fastq_index), max_reads)
    temp_q_count = np.zeros((read_length, len(quality_scores)), dtype=int)
    qual_score_counter = {x: 0 for x in quality_scores}
    quarters = total_records_to_read//4

    i = 0
    wrong_len = 0

    # SeqIO eats up way too much memory for larger fastqs so we're trying to read the file in line by line here
    with open_input(input_file) as fq_in:
        while i < total_records_to_read:

            # We throw away 3 lines and read the 4th, because that's fastq format
            for _ in (0, 1, 2):
                try:
                    fq_in.readline()
                except:
                    break
            line = fq_in.readline()

            """
            This section filters and adjusts the qualities to check. It handles cases of irregular read-lengths as well.
            """
            qualities_to_check = convert_quality_string(line.strip(), qual_offset)

            if len(qualities_to_check) != read_length:
                total_records_to_read += 1
                wrong_len += 1
                if wrong_len % 100 == 0:
                    _LOG.debug(f'So far have detected {wrong_len} reads not matching the mode.')
                continue

            i += 1

            for j in range(read_length):
                # The qualities of each read_position_scores
                quality_bin = take_closest(quality_scores, qualities_to_check[j])
                bin_index = quality_scores.index(quality_bin)
                temp_q_count[j][bin_index] += 1
                qual_score_counter[quality_bin] += 1

            if i % quarters == 0:
                _LOG.info(f'reading data: {(i / total_records_to_read) * 100:.0f}%')

    _LOG.info(f'reading data: 100%')
    _LOG.debug(f'{wrong_len} total reads had a length other than {read_length}')

    avg_std_by_pos = []
    q_count_by_pos = np.asarray(temp_q_count)
    for i in range(read_length):
        this_counts = q_count_by_pos[i]
        expanded_counts = expand_counts(this_counts, quality_scores)
        average_q = np.average(expanded_counts)
        st_d_q = np.std(expanded_counts)
        avg_std_by_pos.append((average_q, st_d_q))

    # Calculates the average error rate
    tot_bases = len(temp_q_count) * len(fastq_index)
    avg_err = 0
    for score in quality_scores:
        error_val = 10. ** (-score / 10.)
        _LOG.info(f"q_score={score}, error value={error_val:e}, count={qual_score_counter[score]}")
        avg_err += error_val * (qual_score_counter[score] / tot_bases)
    _LOG.info(f'Average error rate for dataset: {avg_err}')

    # Generate the sequencing error model with default average error rate
    return avg_std_by_pos, avg_err, read_length
