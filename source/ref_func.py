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


def find_n_regions(reference) -> dict:
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
    # data explanation: my_dat[n_atlas[0][0]:n_atlas[0][1]] = solid block of Ns
    prev_ni = 0
    n_count = 0
    n_atlas = {x: [] for x in reference}
    for chrom in reference:
        sequence = reference[chrom]
        for i in range(len(sequence)):
            if sequence[i] not in ALLOWED_NUCL:
                if n_count == 0:
                    prev_ni = i
                n_count += 1
                if i == len(sequence) - 1:
                    n_atlas[chrom].append((prev_ni, prev_ni + n_count))
            else:
                if n_count > 0:
                    n_atlas[chrom].append((prev_ni, prev_ni + n_count))
                n_count = 0
    non_n_atlas = {x: [] for x in reference}
    for chrom in reference:
        if n_atlas[chrom]:
            # See where N's start (usually at 0)
            if n_atlas[chrom][0][0] > 0:
                # If N's start after 0, then account for that region here
                non_n_atlas[chrom].append((0, n_atlas[chrom][0][0]))
            prev_end = n_atlas[chrom][0][1]
            # skip the first one since we just dealt with it
            for n_region in n_atlas[chrom][1:]:
                non_n_atlas[chrom].append((prev_end, n_region[0]))
                prev_end = n_region[1]
            if prev_end < len(reference[chrom]):
                non_n_atlas[chrom].append((prev_end, len(reference[chrom])))
        else:
            non_n_atlas[chrom].append((0, len(reference[chrom])))
    return non_n_atlas


def model_trinucs(reference, models, safe_zones):
    # Set up the model dicitonary
    trinuc_models = {x: None for x in reference}
    for chrom in reference:
        trinuc_bias = np.zeros(len(reference[chrom]))
        chrom_atlas = safe_zones[chrom]
        non_n = []
        prev_end = 0

        for zone in safe_zones:
            # Start at +1 so we get a trinucleotide to start and end one shy for the same reason
            for i in range(zone[0] + 1, zone[1] - 1):
                trinuc = reference[chrom].seq[i-1:i+2]
                # Let's double check to make sure we didn't pick up a stray N
                if any([i for i in trinuc if i not in ALLOWED_NUCL]):
                    continue
                trinuc_bias[i] = models.mutation_model['trinuc_bias'][trinuc]
        trinuc_models[chrom] = DiscreteDistribution(range(len(reference[chrom])), trinuc_bias)
    return trinuc_models


def process_reference(reference_index, models):
    allowed_areas = find_n_regions(reference_index)
    trinuc_models = model_trinucs(reference_index, models, allowed_areas)
    return allowed_areas, trinuc_models
