"""
Utilities for the sequencing error model
"""

from bisect import bisect_left


__all__ = [
    "bin_scores"
]

DATA_BINS = {}


def bin_scores(bins: list, quality_array: list | str, qual_offset: int):
    """
    Assumes bins list is sorted. Returns the closest value to quality.

    If two numbers are equally close, return the smallest number. Note that in the case of quality scores
    not being binned, the "bin" will end up just being the quality score.

    :param bins: the possible values of the quality scores for the simulation
    :param quality_array: the quality array from data which we will bin.
    :param qual_offset: the
    """
    ret_list = []

    # Convert to array, if string
    if type(quality_array) == str:
        quality_array = convert_to_array(quality_array, qual_offset)

    for score in quality_array:
        if score in DATA_BINS:
            ret_list.append(DATA_BINS[score])
            continue
        pos = bisect_left(bins, score)
        if pos == 0:
            DATA_BINS[score] = bins[0]
            ret_list.append(bins[0])
        elif pos == len(bins):
            DATA_BINS[score] = bins[-1]
            ret_list.append(bins[-1])
        else:
            before = bins[pos - 1]
            after = bins[pos]
            if after - score < score - before:
                DATA_BINS[score] = after
                ret_list.append(after)
            else:
                DATA_BINS[score] = before
                ret_list.append(before)

    return ret_list


def convert_to_array(quality_string: str, off_q: int):
    """
    Converts a quality string from a fastq into the eqivalent array.

    :param quality_string: the quality array in string form
    :param off_q: The offset for the quality str to int conversion
    :return: list of quality scores
    """
    ret_array = []
    for char in quality_string:
        ret_array.append(ord(char) - off_q)

    return ret_array
