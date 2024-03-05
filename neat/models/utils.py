"""
Utilities for the sequencing error model
"""
import numpy as np

from bisect import bisect_left


__all__ = [
    "bin_scores",
    "take_closest"
]

DATA_BINS = {}


def bin_scores(bins: list | np.ndarray, quality_array: list | np.ndarray, qual_offset: int = None):
    """
    Assumes bins list is sorted. Returns the closest value to quality.

    If two numbers are equally close, return the smallest number. Note that in the case of quality scores
    not being binned, the "bin" will end up just being the quality score.

    :param bins: the possible values of the quality scores for the simulation
    :param quality_array: the quality array from data which we will bin.
    :param qual_offset: the quality offset for these quality scores. Only needed if input is a str.
    """
    # Convert to array, if string
    if type(quality_array[0]) == str:
        temp_array = []
        for char in quality_array:
            temp_array.append(ord(char) - qual_offset)
        quality_array = temp_array

    return_array = []
    for score in quality_array:
        if score in DATA_BINS:
            return_array.append(DATA_BINS[score])
        else:
            pos = bisect_left(bins, score)
            if pos == 0:
                DATA_BINS[score] = bins[0]
                return_array.append(bins[0])
            elif pos == len(bins):
                DATA_BINS[score] = bins[-1]
                return_array.append(bins[-1])
            else:
                before = bins[pos - 1]
                after = bins[pos]
                if after - score < score - before:
                    DATA_BINS[score] = after
                    return_array.append(after)
                else:
                    DATA_BINS[score] = before
                    return_array.append(before)

    return return_array

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
