"""
Utilities to generate the sequencing error model
"""

import logging
import numpy as np
import gzip

from Bio import SeqIO

from scipy.stats import mode
from ..models import bin_scores
from ..common import open_input

__all__ = [
    "parse_file"
]

_LOG = logging.getLogger(__name__)


def parse_file(input_file: str, quality_scores: list, max_reads: int, qual_offset: int):
    """
    Parses an individual file for statistics

    :param input_file: The input file to process
    :param quality_scores: A list of potential quality scores
    :param max_reads: Max number of reads to process for this file
    :param qual_offset: The offset for the quality scores when converting to str
    :return:
    """

    _LOG.info(f'file name: {input_file}')

    # SeqIO seems to be having trouble with larger fastq files. I may try just
    # gzip open and read the file. Requires a bit more conversion. Need to convert string to
    # list of scores
    # fastq_index = SeqIO.index(input_file, 'fastq')
    with open_input(input_file) as fastq_in:
        number_records = 0
        read_names = []
        readlens = []
        qualities_to_check = []
        keep_running = True

        while keep_running:
            line1 = fastq_in.readline().strip()
            line2 = fastq_in.readline().strip()
            line3 = fastq_in.readline().strip()
            line4 = fastq_in.readline().strip()
            if not all([line1, line2, line3, line4]):
                keep_running = False
            else:
                number_records += 1
                read_names.append(line1.split(' ')[0].lstrip('@'))
                qual_score = line4.strip()
                qualities_to_check.append(qual_score)
                readlens.append(len(qual_score))

    readlens = np.array(readlens)

    # Using the statistical mode seems like the right approach here. We expect the readlens to be roughly the same.
    readlen_mode = mode(readlens, axis=None, keepdims=False)
    if int(readlen_mode.count) < (0.5 * len(readlens)):
        _LOG.warning("Highly variable read lengths detected. Results may be less than ideal.")
    if int(readlen_mode.count) < 20:
        raise ValueError(f"Dataset is too scarce or inconsistent to make a model. Try a different input.")

    read_length = int(readlen_mode.mode)

    _LOG.debug(f'Read len of {read_length} over {number_records} samples')

    total_records_to_read = min(len(readlens), max_reads)

    # if total_records_to_read > 10e7:
    #     _LOG.warning("Very large dataset. At this time, reading this entire dataset is not feasible. "
    #                  "We will read a sample of the dataset.")
    #     total_records_to_read = int(10e7)

    _LOG.info(f'Reading {total_records_to_read} records out of {len(readlens)}')

    temp_q_count = []
    qual_score_counter = {x: 0 for x in quality_scores}
    quarters = total_records_to_read//4

    _LOG.info("Processing reads...")
    records_skipped = 0
    for i in range(total_records_to_read):
        """
        This section filters and adjusts the qualities to check. It handles cases of irregular read-lengths as well.
        """
        qual_score = qualities_to_check[i]

        if i % quarters == 0 and i != 0:
            _LOG.info(f'reading data: {(i / total_records_to_read) * 100:.0f}%')

        i += 1

        if len(qual_score) != read_length:
            records_skipped += 1
            continue

        quality_bin_list = bin_scores(quality_scores, qual_score, qual_offset)
        for score in quality_bin_list:
            qual_score_counter[score] += 1
        temp_q_count.append(quality_bin_list)

    _LOG.info(f'reading data: 100%')
    _LOG.debug(f'Skipped {records_skipped}/{number_records} records ({1-(records_skipped/number_records):.0%} passed)')

    _LOG.info(f'Building quality-score model for file')
    avg_std_by_pos = []
    q_count_by_pos = np.asarray(temp_q_count).T
    for i in range(read_length):
        average_q = np.average(q_count_by_pos[i])
        st_d_q = np.std(q_count_by_pos[i])
        avg_std_by_pos.append((average_q, st_d_q))

    # Calculates the average error rate
    _LOG.info('Calculating the average error rate for file')
    tot_bases = len(temp_q_count[0]) * len(readlens)
    avg_err = 0
    for score in quality_scores:
        error_val = 10. ** (-score / 10.)
        _LOG.info(f"q_score={score}, error value={error_val:e}, count={qual_score_counter[score]}")
        avg_err += error_val * (qual_score_counter[score] / tot_bases)
    _LOG.info(f'Average error rate for dataset: {avg_err}')

    # Generate the sequencing error model with default average error rate
    return avg_std_by_pos, avg_err, read_length
