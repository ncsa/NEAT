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


def find_habitable_regions(reference):
    habitable_atlas = {x: "" for x in reference}
    for chrom in reference:
        zone_list = [0] * len(reference[chrom])
        sequence = reference[chrom]
        for i in range(len(sequence)):
            if sequence[i] in ALLOWED_NUCL:
                zone_list[i] = 1
        habitable_atlas[chrom] = "".join([str(x) for x in zone_list])
    return habitable_atlas


def model_trinucs(sequence, models, safe_zones):
    # Set up the model dictionary
    trinuc_models = np.zeros(len(sequence))
    # We're going to rewrite this to operate on subslices of the reference
    # To do that we need to know the safe zones of just this region.
    # Start at +1 so we get a trinucleotide to start and end one shy for the same reason
    for i in range(1, len(sequence) - 1):
        trinuc = sequence[i - 1:i + 2]
        # Let's double check to make sure we didn't pick up a stray N or find ourselves
        # outside of the safe zone
        if any([j for j in trinuc if j not in ALLOWED_NUCL])\
                or not int(safe_zones[i]):
            continue
        trinuc_models[i] = models.mutation_model['trinuc_bias'][trinuc]
    trinuc_models = np.array(trinuc_models)
    return DiscreteDistribution(range(len(sequence)), trinuc_models)


nucs = "ACGTN"
dna_string = [random.choice(nucs) for x in range(10000)]
print(dna_string.index('N'))
reference = {'chr1': dna_string}
start = time.time()
for n in range(1000):
    x = find_habitable_regions(reference)
    rand_int = random.randint(0, 9999)
    y = False
    for j in range(100):
        for zone in x['chr1']:
            if zone[0] <= rand_int < zone[1]:
                y = True
                break
print(f"old: {time.time() - start}")

start = time.time()
for n in range(1000):
    x = find_habitable_regions_alt(reference)
    y = False
    for j in range(100):
        if int(x['chr1'][random.randint(0, 9999)]):
            y = True
print(f'new: {time.time() - start}')
