"""
Utilities for modeling the fragment length distribution
"""

import numpy as np
import pysam
from scipy.stats import median_abs_deviation

from .constants import FILTER_MAPQUAL, FILTER_MEDDEV_M

__all__ = [
    "count_frags",
    "filter_lengths"
]


def count_frags(file: str) -> list:
    """
    Takes a sam or bam file_list input and creates a list of the number of reads that are paired,
    first in the pair, confidently mapped and whose pair is mapped to the same reference
    :param file: A sam input file_list
    :return: A list of the tlens from the bam/sam file_list
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
