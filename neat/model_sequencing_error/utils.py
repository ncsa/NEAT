"""
Utilities to generate the sequencing error model
"""

import logging
import numpy as np

from Bio import SeqIO

from pathlib import Path
from bisect import bisect_left
from scipy.stats import mode

__all__ = [
    "take_closest",
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


def parse_file(input_file: str, quality_scores: list, max_reads: int):
    """
    Parses an individual file for statistics

    :param input_file: The input file to process
    :param quality_scores: A list of potential quality scores
    :param max_reads: Max number of reads to process for this file
    :return:
    """

    _LOG.info(f'reading {input_file}')

    fastq_index = SeqIO.index(input_file, 'fastq')
    read_names = list(fastq_index)

    readlens = []

    counter = 0
    for read_name in fastq_index:
        read = fastq_index[read_name]
        if read.letter_annotations['phred_quality']:
            readlens.append(len(read))
            counter += 1
            if counter > 1000:
                # takes too long and uses too much memory to read all of them, so let's just get a sample.
                break

    readlens = np.array(readlens)

    # Using the statistical mode seems like the right approach here. We expect the readlens to be roughly the same.
    readlen_mode = mode(readlens, axis=None, keepdims=False)
    if readlen_mode.count < (0.5 * len(readlens)):
        _LOG.warning("Highly variable read lengths detected. Results may be less than ideal.")
    if readlen_mode.count < 20:
        raise ValueError(f"Dataset is too scarce or inconsistent to make a model. Try a different input.")

    read_length = readlen_mode.mode

    _LOG.debug(f'Read len of {read_length}, with {readlen_mode.count} out of {len(fastq_index)}')

    # In order to account for scarce data, we may try to set the minimum count at 1.
    # For large datasets, this will have minimal impact, but it will ensure for small datasets
    # that we don't end up with probabilities of scores being 0. To do this, just change np.zeros to np.ones
    temp_q_count = np.zeros((read_length, len(quality_scores)), dtype=int)
    qual_score_counter = {x: 0 for x in quality_scores}

    total_records_to_read = min(len(fastq_index), max_reads)
    quarters = total_records_to_read//4

    rng = np.random.default_rng()

    i = 0
    wrong_len = 0
    while i < total_records_to_read:
        """
        This section filters and adjusts the qualities to check. It handles cases of irregular read-lengths as well.
        """
        # Get the ith key.
        read = fastq_index[read_names[i]]
        qualities_to_check = read.letter_annotations['phred_quality']
        i += 1
        if len(qualities_to_check) != read_length:
            total_records_to_read += 1
            wrong_len += 1
            if wrong_len % 100 == 0:
                _LOG.debug(f'So far have detected {wrong_len} reads not matching the mode.')
            continue

        for j in range(read_length):
            # The qualities of each read_position_scores
            quality_bin = take_closest(quality_scores, qualities_to_check[j])
            bin_index = quality_scores.index(quality_bin)
            temp_q_count[j][bin_index] += 1
            qual_score_counter[quality_bin] += 1

        if i % quarters == 0:
            _LOG.info(f'reading data: {(i / total_records_to_read) * 100:.0f}%')

    _LOG.info(f'reading data: 100%')

    quality_score_probabilities = np.zeros((read_length, len(quality_scores)), dtype=float)

    for i in range(read_length):
        # TODO Add some interpolation for missing scores.
        read_position_scores = temp_q_count[i]
        total = sum(read_position_scores)
        # If total is zero, we had no data at that position, for some reason. Assume a uniform chance of any score.
        if total == 0:
            # No reads at that position for the entire dataset, assume uniform probability
            quality_score_probabilities[i] = np.full(
                len(quality_scores), 1/len(quality_scores), dtype=float
            )
            continue
        if total < 20:
            # fill it out with random scores
            num = 20 - total
            # We'll make sure we have at least 1 'score' at each position, to give us interesting results
            # with the simulation, while keeping the distribution roughly the same.
            temp_read_pos = read_position_scores + 1
            temp_probs = temp_read_pos/sum(temp_read_pos)
            scores = rng.choice(quality_scores, size=num, p=temp_probs)
            read_position_scores = np.concatenate((read_position_scores, scores))
            total = sum(read_position_scores)

        quality_score_probabilities[i] = read_position_scores/total

    # Calculates the average error rate
    tot_bases = float(sum(qual_score_counter.values()))
    avg_err = 0
    for k in sorted(qual_score_counter.keys()):
        error_val = 10. ** (-k / 10.)
        _LOG.info(f"q_score={k}, error value={error_val:e}, count={qual_score_counter[k]}")
        avg_err += error_val * (qual_score_counter[k] / tot_bases)
    _LOG.info(f'Average error rate for dataset: {avg_err}')

    # Generate the sequencing error model with default average error rate
    return quality_score_probabilities, avg_err, read_length
