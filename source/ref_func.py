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


def find_habitable_regions(reference) -> dict:
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


def model_trinucs(reference, models, safe_zones):
    # Set up the model dicitonary
    trinuc_models = {x: None for x in reference}
    for chrom in reference:
        trinuc_bias = np.zeros(len(reference[chrom]))
        chrom_atlas = safe_zones[chrom]

        for zone in chrom_atlas:
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
    allowed_areas = find_habitable_regions(reference_index)
    trinuc_models = model_trinucs(reference_index, models, allowed_areas)
    return allowed_areas, trinuc_models
