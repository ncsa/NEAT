import gzip
import logging
import pathlib
import random
import time
import numpy as np

from Bio.Seq import MutableSeq
from Bio.Seq import Seq

from source.constants_and_defaults import ALLOWED_NUCL
from source.error_handling import premature_exit, log_mssg
from source.probability import DiscreteDistribution


def find_habitable_regions_alt(reference) -> dict:
    """
    Finds N regions in the sequence
    :param reference: Biopython SeqRecord object containing the sequence to scan.
    :return: atlas of n locations in dicitonary form

    >>> my_ref = {'chr1': Seq("NNNNNAAACCCTTTNAAAACCCCNNNNN")}
    >>> find_n_regions(my_ref)
    {'chr1': [(5, 14), (15, 23)]}

    >>> my_ref = {'chr1': Seq("ACATGANNNNNAAACCCTTTNAAAACCCCNNNNN")}
    >>> find_n_regions(my_ref)
    {'chr1': [(0, 6), (11, 20), (21, 29)]}

    >>> my_ref = {'chr1': Seq("ACATGAAAACCCTTTAAAACCCC")}
    >>> find_n_regions(my_ref)
    {'chr1': [(0, 23)]}
    """
    # TODO instead, this may need to be per chromosome OR written to a bed file, for memory reasons, though I think
    #  the final product is not too big
    # data explanation: my_dat[n_atlas[0][0]:n_atlas[0][1]] = solid block of Ns
    non_n_atlas = {x: [] for x in reference}
    for chrom in reference:
        sequence = reference[chrom]
        prev_non_ni = 0
        non_n_count = 0
        for i in range(len(sequence)):
            if sequence[i] in ALLOWED_NUCL:
                if non_n_count == 0:
                    prev_non_ni = i
                non_n_count += 1
                if i == len(sequence) - 1:
                    non_n_atlas[chrom].append((prev_non_ni, prev_non_ni + non_n_count))
            else:
                if non_n_count > 0:
                    non_n_atlas[chrom].append((prev_non_ni, prev_non_ni + non_n_count))
                non_n_count = 0
    return non_n_atlas


