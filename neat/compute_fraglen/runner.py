"""
Creates a model of fragment lengths in a dataset
"""

import gzip
import pickle
import numpy as np
import pysam
import logging

from pathlib import Path

from ..models import FragmentLengthModel
from .constants import FILTER_MAPQUAL
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


def compute_fraglen_runner(file: str | Path, output: str | Path):
    """
    Main function takes 2 arguments:

    :param file: a path to a sam or bam file input.
    :param output: the string prefix of the output
    """
    input_file = file
    output_prefix = output
    output = Path(output_prefix + '.pickle.gz')
    validate_output_path(output)

    all_tlens = count_frags(input_file)
    if not all_tlens:
        raise ValueError("No valid template lengths in sam file.")

    st_dev = float(np.std(all_tlens))
    mean = float(np.mean(all_tlens))
    max_tlen = max(all_tlens)
    min_tlen = min(all_tlens)

    model = FragmentLengthModel(st_dev, mean, max_tlen, min_tlen)
    print('\nSaving model...')
    with gzip.open(output) as outfile:
        pickle.dump(model, outfile)
