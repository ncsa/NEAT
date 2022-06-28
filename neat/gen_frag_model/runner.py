"""
Creates a model of fragment lengths in a dataset
"""

import gzip
import pickle
import numpy as np
import pysam
import logging
from scipy.stats import median_abs_deviation

from pathlib import Path

from ..models import FragmentLengthModel
from .constants import FILTER_MAPQUAL, FILTER_MEDDEV_M
from ..common import validate_output_path

__all__ = [
    "compute_fraglen_runner"
]

_LOG = logging.getLogger(__name__)


def count_frags(file: str) -> list:
    """
    Takes a sam or bam file input and creates a list of the number of reads that are paired,
    first in the pair, confidently mapped and whose pair is mapped to the same reference
    :param file: A sam input file
    :return: A list of the tlens from the bam/sam file
    """
    count_list = []

    with pysam.AlignmentFile(file) as file_to_parse:
        for item in file_to_parse:
            # new values based on pysam
            sam_flag = item.flag
            my_ref = item.reference_id
            map_qual = item.mapping_quality
            mate_ref = item.next_reference_id
            my_tlen = abs(item.template_length)

            # if read is paired, and is first in pair, and is confidently mapped...
            if sam_flag & 1 and sam_flag & 64 and map_qual > FILTER_MAPQUAL:
                # and mate is mapped to same reference
                if mate_ref == '=' or mate_ref == my_ref:
                    count_list.append(my_tlen)

    return count_list


def filter_lengths(datalist: list, min_reads: int) -> list:
    """
    Filters the datalist to remove outliers.

    TODO this function might be useful in the simulation as well

    :param datalist: The list to filter
    :param min_reads: filter reads with less than this many reads. If 0, filter_lengths returns the original datalist.
    :return: The filtered list
    """
    # If min_reads set to 0, then skip filtering.
    if not min_reads:
        return datalist

    med = np.median(datalist)
    mad = median_abs_deviation(datalist)
    output_list = []
    for item in (list(set(datalist))):
        if 0 < item <= med + FILTER_MEDDEV_M * mad:
            data_count = datalist.count(item)
            if data_count >= min_reads:
                output_list.extend([item] * data_count)

    return output_list


def compute_fraglen_runner(file: str | Path, filter_minreads: int, output: str | Path):
    """
    Main function takes 2 arguments:

    :param file: a path to a sam or bam file input.
    :param filter_minreads: minimum number of reads needed to count a fragment length.
                            If 0, then filtering will be skipped
    :param output: the string prefix of the output
    """
    _LOG.info("Generating fragment length model")

    input_file = file
    output_prefix = output
    output = Path(output_prefix + '.pickle.gz')
    validate_output_path(output)

    all_tlens = count_frags(input_file)
    if not all_tlens:
        raise ValueError("No valid template lengths in sam file.")

    filtered_lengths = filter_lengths(all_tlens, filter_minreads)

    if not filtered_lengths:
        raise ValueError("No data passed the filter, nothing to calculate. Try adjusting the filter settings.")

    st_dev = float(np.std(filtered_lengths))
    mean = float(np.mean(filtered_lengths))
    max_tlen = max(filtered_lengths)
    min_tlen = min(filtered_lengths)

    model = FragmentLengthModel(st_dev, mean, max_tlen, min_tlen)
    print('\nSaving model...')
    with gzip.open(output, 'w+') as outfile:
        pickle.dump(model, outfile)
