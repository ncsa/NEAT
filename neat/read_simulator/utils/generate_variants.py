"""
This task generates the variants that will become the VCF or the basis for downstream files. This is the heart
of the NEAT generate_reads task and will always be run in a NEAT generate_reads call.

"""

import logging
import time
import bisect
import numpy as np
import re

from copy import deepcopy
from Bio import SeqRecord
from numpy.random import Generator
from pathlib import Path

from ...models import MutationModel, InsertionModel, DeletionModel, SnvModel
from ...variants import Insertion, Deletion, SingleNucleotideVariant, ContigVariants
from ..utils import Options
from ...common import ALLOWED_NUCL, pick_ploids, NUC_IND, DINUC_IND, open_input, open_output

_LOG = logging.getLogger(__name__)

__all__ = [
    "generate_variants"
]


def find_random_non_n(rng: Generator,
                      safe_zones: np.ndarray):

    # Normalize the n-gap array to create a probability array
    safe_zones = safe_zones / sum(safe_zones)
    positions = np.arange(len(safe_zones))
    return int(rng.choice(positions, p=safe_zones))


def map_non_n_regions(sequence) -> np.ndarray:
    """
    We know our endpoints are both allowed characters, so we just need to make sure
    that the sequence in question isn't too N-heavy
    :param sequence: the sequence to analyze
    :return: None if N-concentration is over 10%, the non-n map otherwise
    """
    base_check = np.ones(len(sequence), dtype=int)
    for gap in re.finditer("N+", str(sequence)):
        base_check[gap.start(): gap.end()] = 0

    # Less than 90% non-N's we'll skip
    if np.average(base_check) <= 0.90:
        return np.empty(0)

    return base_check


def generate_variants(reference: SeqRecord,
                      mutation_rate_regions: list[tuple[int, int]],
                      existing_variants: ContigVariants,
                      mutation_model: MutationModel,
                      options: Options,
                      max_qual_score: int):
    """
    This function will generate variants to add to the dataset, by writing them to the input temp vcf file.

    TODO: need to add cancer logic to this section

    :param reference: The reference to generate variants off of. This should be a SeqRecord object
    :param mutation_rate_regions: Genome segments with associated mutation rates
    :param existing_variants: Any input variants or overlaps to pick up from the previous contig
    :param mutation_model: The mutation model for the dataset
    :param options: The options for the run
    :param max_qual_score: the maximum quality score for the run.
    """
    return_variants = existing_variants
    # Step 1: Create a VCF of mutations

    # pase out the mutation rates
    mutation_rates = np.array([x[2] for x in mutation_rate_regions])

    # Trying to use a random window to keep memory under control. May need to adjust this number.
    max_window_size = 1000
    n_added = 0

    start = time.time()

    """
    Decide how many and where any random mutations will happen in this contig
    """
    total = len(reference)
    _LOG.debug(f"Contig length: {total}")
    factors = []
    for i in range(len(mutation_rate_regions)):
        region_len = abs(mutation_rate_regions[i][1] - mutation_rate_regions[i][0])
        # Weigh each region by its total length
        factors.append(region_len * mutation_rates[i])
        total -= region_len
    # Whatever is left is weighted by the average rate
    factors.append(total * mutation_model.avg_mut_rate)
    overall_mutation_average = sum(factors)/len(reference)
    average_number_of_mutations = round(len(reference) * overall_mutation_average)
    max_mutations = round(len(reference) * 0.3)

    # This number will serve as a counter for our loops. Default is 1 mutation per chromosome.
    min_mutations = options.min_mutations

    _LOG.debug(f"Average number of mutations for this contig: {average_number_of_mutations}")
    # Pick a random number from a poisson distribution. We want at least min_mutations and at most max mutations.
    how_many_mutations = min(max(options.rng.poisson(average_number_of_mutations), min_mutations), max_mutations)

    _LOG.debug(f'Planning to add {how_many_mutations} mutations. The final number may be less.')

    while how_many_mutations > 0:
        # Pick a region based on the mutation rates
        # (default is one rate for the whole chromosome, so this will be trivial in that case
        # for this selection, we'll normalize the mutation rates
        probability_rates = mutation_rates / sum(mutation_rates)
        mut_region = options.rng.choice(a=mutation_rate_regions, p=probability_rates)
        # Pick a random starting place. Randint is inclusive of endpoints, so we subtract 1
        window_start = options.rng.integers(mut_region[0], mut_region[1] - 1, dtype=int)
        found = False
        if reference[window_start] not in ALLOWED_NUCL:
            # pick a random location to the right
            plus = options.rng.integers(window_start + 1, mut_region[1] - 1, dtype=int)
            if reference[plus] in ALLOWED_NUCL:
                found = True
                window_start = plus
            else:
                # If that didn't work pick a random location to the left
                if window_start - 1 > mut_region[0]:
                    minus = options.rng.integers(mut_region[0], window_start - 1, dtype=int)
                    if reference[minus] in ALLOWED_NUCL:
                        found = True
                        window_start = minus
        else:
            found = True

        # If we couldn't find a spot, try again
        if not found:
            continue

        # at this point the location is found. Now we need a second endpoint. Grab 1000 bases or to the end.
        end_point = options.rng.integers(window_start,
                                         min(mut_region[1], window_start + max_window_size) - 1,
                                         dtype=int)
        if reference[end_point] not in ALLOWED_NUCL:
            # Didn't find it to the right, look to the left
            end_point = options.rng.integers(max(window_start - max_window_size, mut_region[0]),
                                             window_start,
                                             dtype=int)
            if reference[end_point] not in ALLOWED_NUCL:
                # No suitable end_point, so we try again
                continue

        # Sorting assures that wherever we found the end point, the coordinates will be in the correct order for slicing
        mutation_slice = sorted([window_start, end_point])
        slice_distance = mutation_slice[1] - mutation_slice[0]
        # How many variants to add in this slice. currently set to 1, because if we got this far we have at least
        # 1 mutation to add and if we set this to 0, it might spin all day before it finds a suitable location.
        variants_to_add_in_slice = max(int((slice_distance/len(reference)) * how_many_mutations), 1)

        subsequence = reference[mutation_slice[0]: mutation_slice[1]].seq.upper()
        # In the case where the end points are equal or the segment is too small, just skip
        if len(subsequence) <= 10:
            continue

        # We'll find the n-gaps and for the slice.
        n_gaps = map_non_n_regions(subsequence)

        # If the sequence has too many N's, we'll skip it
        if not any(n_gaps):
            continue

        # Now we model trinuclotide bias, for SNVs in the slice
        mutation_model.map_local_trinuc_bias(subsequence, n_gaps)

        # Begin random mutations for this slice
        # Note that any new variant types will need code in this area to handle the functions.
        debug = 0
        while variants_to_add_in_slice > 0:
            # We decrement now because we don't want to get stuck in a never ending loop
            variants_to_add_in_slice -= 1
            how_many_mutations -= 1

            # Now figure out the type of random mutation to insert
            variant_type = mutation_model.get_mutation_type()

            # Case 1: indel
            if variant_type == Insertion or variant_type == Deletion:
                # Because the slice is a substring, we add the start point onto the relative location to get the
                # location relative to the reference.
                position = find_random_non_n(options.rng, n_gaps)  # position in slice
                location = position + mutation_slice[0]  # location relative to reference
                if variant_type == Insertion:
                    temp_variant = mutation_model.generate_insertion(location, subsequence[position])
                else:
                    temp_variant = mutation_model.generate_deletion(location)

            # Case 2: SNV
            elif variant_type == SingleNucleotideVariant:
                # We'll sample for the location within this slice
                # It's a relative location, so we add the start point of the subsequence to that.
                position = mutation_model.sample_trinucs()  # position in slice
                location = position + mutation_slice[0]  # location relative to reference
                if location == 0:
                    continue
                trinuc = reference[location: location+3].seq.upper()
                temp_variant = mutation_model.generate_snv(trinuc, location)

            else:
                raise ValueError(f"Attempting to create an unsupported variant: {variant_type}")

            # pick which ploid is mutated
            temp_variant.genotype = pick_ploids(options.ploidy, mutation_model.homozygous_freq, 1, options.rng)

            # There shouldn't be a ton of overlapping variants, but this is to handle those.
            if location in return_variants:
                """
                If the location already exists, then we'll need to force it to pick a ploid 
                that currently doesn't have a variant. This overrides the default genotype
                variable created above, but it shouldn't happen very often.
                """
                if return_variants.find_dups(temp_variant):
                    # This compiles all the variants at this location, giving a 1 for every ploid that has a variant.
                    composite_genotype = return_variants.compile_genotypes_for_location(location)
                    if 0 not in composite_genotype:
                        # Here's a counter to make sure we're not getting stuck on a single location
                        debug += 1
                        if debug > 1000000:
                            raise RuntimeError("Check this if, as it may be causing an infinite loop.")
                        # No suitable place to put this, so we skip.
                        continue
                    # This sets up a probability array with weights 1 for open spots (x==0) and 0 elsewhere
                    probs = np.array([1 if x == 0 else 0 for x in composite_genotype])
                    probs = probs / sum(probs)
                    # Pick an index of a position to mutate based on the probabilities, which are uniform for 0s left
                    # in the composite genotype
                    ploid = options.rng.choice(list(range(len(composite_genotype))), p=probs)
                    genotype = np.zeros(options.ploidy)
                    genotype[ploid] = 1
                    temp_variant.genotype = genotype

            # Make sure this new variant doesn't overlap an existing insertion or deletion
            in_deletion = return_variants.check_if_del(temp_variant)
            in_insertion = return_variants.check_if_ins(temp_variant)

            if in_deletion:
                if type(temp_variant) == Insertion:
                    temp_variant.position1 = in_deletion.position1
                else:
                    continue

            if in_insertion:
                if type(temp_variant) == Deletion:
                    # this is a delins, so we set the location to the same as the deletion
                    temp_variant.position1 = in_insertion.position1
                else:
                    continue

            temp_variant.qual_score = max_qual_score
            return_variants.add_variant(temp_variant)
            # The count will tell us how many we actually added v how many we were trying to add
            n_added += 1

    _LOG.info(f'Finished generating random mutations in {(time.time() - start)/60:.2f} minutes')
    _LOG.info(f'Added {n_added} mutations to {reference.id}')

    return return_variants


