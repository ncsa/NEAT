"""
Utilities to generate the sequencing error model
"""

import logging
import numpy as np

from Bio import SeqIO

from pathlib import Path
import pickle
from bisect import bisect_left
from scipy.stats import mode

__all__ = [
    "bin_scores",
    "parse_file"
]

_LOG = logging.getLogger(__name__)


def bin_scores(bins, quality_array):
    """
    Assumes bins list is sorted. Returns the closest value to quality.

    If two numbers are equally close, return the smallest number. Note that in the case of quality scores
    not being binned, the "bin" will end up just being the quality score.

    :bins list: the possible values of the quality scores for the simulation
    :quality_array list: the quality array from data which we will bin.
    """
    ret_list = []

    for score in quality_array:
        pos = bisect_left(bins, score)
        if pos == 0:
            ret_list.append(bins[0])
        elif pos == len(bins):
            ret_list.append(bins[-1])
        else:
            before = bins[pos - 1]
            after = bins[pos]
            if after - score < score - before:
                ret_list.append(after)
            else:
                ret_list.append(before)

    return ret_list


def filter_reads(my_reads, length):
    """
    Filters a list of reads down to speed up processing

    :param my_reads: The fastq index of reads which we will filter
    :param length: The length that we will use as the filter parameter
    """
    records_to_return = []

    records_skipped = 0
    for record in my_reads:
        qual_score_len = len(my_reads[record].letter_annotations['phred_quality'])
        if qual_score_len < length:
            records_skipped += 1
            continue
        elif qual_score_len == length:
            records_to_return.append(record)
        else:
            # We could try truncating these down. Or we could skip them. Not sure which is better.
            records_skipped += 1
            continue

    _LOG.info(f'{records_skipped} out of {len(my_reads)} were skipped for having a weird read length')

    return records_to_return


def parse_file(input_file: str, quality_scores: list, max_reads: int):
    """
    Parses an individual file for statistics

    :param input_file: The input file to process
    :param quality_scores: A list of potential quality scores
    :param max_reads: Max number of reads to process for this file
    :return:
    """

    _LOG.info(f'file name: {input_file}')

    fastq_index = SeqIO.index(input_file, 'fastq')
    number_records = len(fastq_index)
    read_names = list(fastq_index)

    readlens = []

    counter = 0
    for read_name in fastq_index:
        read = fastq_index[read_name]
        if read.letter_annotations['phred_quality']:
            readlens.append(len(read))
            counter += 1
            if counter >= number_records//100:
                # takes too long and uses too much memory to read all of them, so let's just get a 1% sample.
                break

    readlens = np.array(readlens)

    # Using the statistical mode seems like the right approach here. We expect the readlens to be roughly the same.
    readlen_mode = mode(readlens, axis=None, keepdims=False)
    if int(readlen_mode.count) < (0.5 * len(readlens)):
        _LOG.warning("Highly variable read lengths detected. Results may be less than ideal.")
    if int(readlen_mode.count) < 20:
        raise ValueError(f"Dataset is too scarce or inconsistent to make a model. Try a different input.")

    read_length = int(readlen_mode.mode)

    _LOG.debug(f'Read len of {read_length} over {counter} samples')

    total_records_to_read = min(len(fastq_index), max_reads)
    temp_q_count = []
    qual_score_counter = {x: 0 for x in quality_scores}
    quarters = total_records_to_read//4

    i = 0
    records_skipped = 0
    while i < total_records_to_read:
        """
        This section filters and adjusts the qualities to check. It handles cases of irregular read-lengths as well.
        """
        # Get the ith key.
        read = fastq_index[read_names[i]]
        qualities_to_check = read.letter_annotations['phred_quality']

        if len(qualities_to_check) != read_length:
            records_skipped += 1
            continue

        quality_bin_list = bin_scores(quality_scores, qualities_to_check)
        for score in quality_bin_list:
            qual_score_counter[score] += 1
        temp_q_count.append(quality_bin_list)

        i += 1

        if i % quarters == 0:
            _LOG.info(f'reading data: {(i / total_records_to_read) * 100:.0f}%')

    _LOG.info(f'reading data: 100%')
    _LOG.debug(f'Skipped {records_skipped}/{number_records} records ({1-(records_skipped/number_records):.0%}% passed)')

    _LOG.info(f'Building quality-score model for file')
    avg_std_by_pos = []
    q_count_by_pos = np.asarray(temp_q_count).T
    for i in range(read_length):
        average_q = np.average(q_count_by_pos[i])
        st_d_q = np.std(q_count_by_pos[i])
        avg_std_by_pos.append((average_q, st_d_q))

    # Calculates the average error rate
    _LOG.info('Calculating the average error rate for file')
    tot_bases = len(temp_q_count[0]) * len(fastq_index)
    avg_err = 0
    for score in quality_scores:
        error_val = 10. ** (-score / 10.)
        _LOG.info(f"q_score={score}, error value={error_val:e}, count={qual_score_counter[score]}")
        avg_err += error_val * (qual_score_counter[score] / tot_bases)
    _LOG.info(f'Average error rate for dataset: {avg_err}')

    # Generate the sequencing error model with default average error rate
    return avg_std_by_pos, avg_err, read_length
