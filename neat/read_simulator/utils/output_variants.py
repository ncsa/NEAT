"""
This task generates the variants that will become the VCF or the basis for downstream files. This is the heart
of the NEAT generate_reads task and will always be run in a NEAT generate_reads call.
"""

import logging
import numpy as np
import time
import bisect

from Bio import SeqRecord
from numpy.random import Generator

import neat.models as models
from .options import Options
from .ploid_functions import pick_ploids
from ...common import ALLOWED_NUCL, open_input, open_output

_LOG = logging.getLogger(__name__)


def add_variant_to_fasta(position: int, ref: str, alt: str, reference_to_alter: list):
    ref_len = len(ref)
    # Apply mutation
    reference_to_alter[position:position+ref_len] = alt
    return reference_to_alter


def check_if_deleted(all_dels, location):

    for deletion in all_dels:
        # deletion[0] = location of deletion
        # deletion[1] = length of the deletion
        if deletion[0] < location < deletion[0] + deletion[1]:
            return True

    return False


def find_random_non_n(rng: Generator, mutation_slice, safe_zones, max_attempts):
    for _ in range(max_attempts):
        potential_location = rng.integers(0, len(safe_zones), dtype=int)
        if not safe_zones[potential_location]:
            continue
        else:
            return potential_location + mutation_slice[0]

    return False


def map_non_n_regions(sequence):
    """
    We know our endpoints are both allowed characters, so we just need to make sure
    that the sequence in question isn't too N-heavy
    :param sequence: the sequence to analyze
    :return: False if N-concentration is over 25%, true otherwise
    """
    base_check = [0] * len(sequence)
    for i in range(len(sequence)):
        if sequence[i] in ALLOWED_NUCL:
            base_check[i] = 1

    # Greater than 10% N's we'll skip
    if sum(base_check)/len(sequence) <= 0.90:
        return None

    return base_check


def parse_mutation_rate_dict(avg_rate: float,
                             reference: SeqRecord,
                             mutation_rate_map: dict[dict[int, int], dict[float]] = None)\
                             -> (list[int, int], list[float]):
    """
    This parses the mutation rate dict and fills in any gaps, so it can be more easily cycled through
    later.

    example mutation_rate_map {'H1N1_HA': {'regions': [(22, 500), (510, 750)],
                                           'rates': [0.001, , 0.003]}
                               'H1N1_PB': {...}, ...}
    The input to this function is the dict for a single chromosome.

    :param avg_rate: The average mutation rate desired for this dataset
    :param reference: The SeqRecord for this contig
    :param mutation_rate_map: A dictionary based on the input from the bed file. If there was no input,
                              then this will set a uniform mutation rate across the chromosome.
    :return: A tuple with the list of regions (tuples with start, end indexes)
             and a list of the corresponing mutation rates
    """

    region_list = []
    mutation_list = []
    start = 0
    if mutation_rate_map:
        for region in mutation_rate_map:
            if region[0] > start:
                region_list.append((start, region[0]))
                mutation_list.append(avg_rate)
                start = region[1]
                region_list.append((region[0], region[1]))
                mutation_list.append(region[2])
            elif region[0] == start:
                region_list.append((region[0], region[1]))
                mutation_list.append(region[2])
                start = region[1]
    else:
        region_list.append((0, len(reference)))
        mutation_list.append(avg_rate)
    if region_list[-1][1] != len(reference):
        region_list.append((start, len(reference)))
        mutation_list.append(avg_rate)
    return region_list, mutation_list


def generate_variants(reference: SeqRecord,
                      tmp_vcf_fn: str,
                      trinucleotide_map: list,
                      target_regions: dict,
                      discard_regions: dict,
                      mutation_rate_regions: dict[dict[int, int], dict[float]],
                      output_variants: dict,
                      error_model: models.SequencingErrorModel,
                      mut_model: models.MutationModel,
                      options: Options):
    """
    This function will generate variants to add to the dataset. In the event that the user only wants a vcf, then
    this will translate the vcf.

    TODO: rewrite variants so they have a consistent format, in order to facilitate structural variants
     it makes sense to me to do these in the temporary vcf

    One representation of variants:
    DUP: pos, len
    INV: pos, len
    DEL: pos, len
    INS: pos1, len, pos2,
    TRA: pos1, len1, pos2, len2

    Here is the same thing done in a consistent format:
    DUP: pos1, len1, pos2=pos1+len1, len2=NULL
    INV: pos1, len1, pos2=NULL, len2=-len1
    DEL: pos1, len1, pos2=NULL, len2=NULL
    INS: pos1, len1, pos2, len2=NULL
    TRA: pos1, len1, pos2, len2

    TODO: need to add cancer logic to this section

    :param reference: A record from a fasta file
    :param tmp_vcf_fn: The filename where to write the temp vcf for this chromosome
    :param trinucleotide_map: The trinucleotide map of this chromosome
    :param target_regions: User input regions to target
    :param discard_regions: User input regions to discard
    :param mutation_rate_regions: User input regions with mutation rates, parsed into a dictionary
    :param output_variants: Variants input by the user that must be included in the output
    :param error_model: The sequencing error model derived from real data
    :param mut_model: The mutation model derived from real data
    :param options: User options for this run
    :return:
        tmp_vcf_fn: The NEAT generated mutations in vcf format
        mutated_references: The mutated fasta or fastas for downstream use
    """

    # Step 1: Create a VCF of mutations
    _LOG.info(f'Generating chromosome mutation map.')
    mutation_map, mutation_rates = parse_mutation_rate_dict(mut_model.avg_mut_rate, reference, mutation_rate_regions)

    # Trying to use a random window to keep memory under control. May need to adjust this number.
    max_window_size = 1000
    n_added = 0
    random_mutation_quality_score = max(error_model.quality_scores)
    random_mutation_filter = "PASS"
    """
    Checking current mutations from input vcf. If there were no variants, output_variants will be empty
    """
    output_variants_locations = []
    if output_variants:
        """
        For reference, the columns in output_variants, and their indices:
        key: POS
        values[0]: ID [0], REF [1], ALT [2], QUAL [3], FILTER [4], INFO [5], 
                FORMAT [6, optional], SAMPLE1 [7, optional], SAMPLE2 [8, optional]
        values[1]: genotype (if present in original vcf), None if not present
        values[2]: genotype tumor (if present in original vcf), None if not present and a cancer sample
        """
        # this should give is just the positions of the inserts. This will be used to keep track of sort order.
        output_variants_locations = sorted(list(output_variants.keys()))

    _LOG.info(f'Adding random mutations for {reference.id}')
    start = time.time()

    """
    Decide how many and where any random mutations will happen in this contig
    """
    total = len(reference)
    factors = []
    for region in mutation_rate_regions:
        region_len = abs(region[1] - region[0])
        # Weigh each region by its total length
        factors.append(region_len * region[2])
        total -= region_len
    # Whatever is left is weighted by the average rate
    factors.append(total * mut_model.avg_mut_rate)
    overall_mutation_average = sum(factors)/len(reference)
    average_number_of_mutations = round(len(reference) * overall_mutation_average)
    max_mutations = round(len(reference) * 0.3)

    # This number will serve as a counter for our loops. Default is 1 mutation per chromosome.
    min_mutations = options.min_mutations

    # Pick a random number from a poisson distribution. We want at least min_mutations and at most max mutations.
    how_many_mutations = min(max(options.rng.poisson(average_number_of_mutations), min_mutations), max_mutations)

    _LOG.info(f'Planning to add {how_many_mutations} mutations. The final number may be less.')

    # We may need to skip locations if they were deleted
    all_dels = []

    while how_many_mutations > 0:
        # Pick a region based on the mutation rates
        # (default is one rate for the whole chromosome, so this will be trivial in that case
        mut_region = options.rng.choice(a=mutation_map, p=mutation_rates)
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

        subsequence = reference[mutation_slice[0]: mutation_slice[1]]
        # In the case where the end points are equal, just skip
        if len(subsequence) == 0:
            continue

        non_n_regions = map_non_n_regions(subsequence)

        # If the sequence has too many N's, we'll skip it
        if not non_n_regions:
            continue

        # Begin random mutations for this slice
        while variants_to_add_in_slice > 0:
            # We decrement now because we don't want to get stuck in a never ending loop
            variants_to_add_in_slice -= 1
            how_many_mutations -= 1

            # Now figure out the type of random mutation to insert
            variant_type = mut_model.get_mutation_type()
            # Case 1: indel
            if variant_type == models.Insertion or variant_type == models.Deletion:
                # First pick a location. This function ensures that we do not pick an N as our starting place
                # Note that find_random_non_n helpfully adds the slice start to the location.
                location = find_random_non_n(options.rng, mutation_slice, non_n_regions, 5)

                if check_if_deleted(all_dels, location):
                    continue  # No increments, no attempts, just try again.

                if variant_type == models.Insertion:
                    length = mut_model.get_insertion_length()
                    # Try to find a location to insert. Give it ten tries, then give up.
                    ref = reference[location]
                    # Check if a p= parameter is needed here
                    insertion = ''.join(options.rng.choice(ALLOWED_NUCL, size=length))
                    alt = ref + insertion
                else:
                    length = mut_model.get_deletion_length()
                    # Plus one so we make sure to grab the first base too
                    ref = reference[location: location+length+1].seq
                    alt = reference[location]
                    all_dels.append((location, length))

            # Case 2: SNP
            else:
                # We'll sample for the location within this slice
                # It's a relative location, so we add the start point of the subsequence to that.
                location = mut_model.sample_trinucs(trinucleotide_map[mutation_slice[0]: mutation_slice[1]]) \
                           + mutation_slice[0]

                if check_if_deleted(all_dels, location):
                    continue  # No increments, no attempts, just try again.

                ref = reference[location]
                trinuc = reference[location - 1: location + 2]
                transition_values = ALLOWED_NUCL
                # TODO figure out how to fix this next line. Pretty sure this won't work as written
                transition_probs = list(mut_model.trinuc_trans_matrices[trinuc.seq])
                # We want this to be an actual variant, so we'll try a few times
                alt = options.rng.choice(transition_values, p=transition_probs)
                # Max 10 tries
                j = 10
                while alt == ref and j > 0:
                    alt = options.rng.choice(transition_values, p=transition_probs)
                    j -= 1

                # We tried and failed, so we give up.
                if alt == ref:
                    continue

            # Now we have a location, an alt and a ref, so we figure out which ploid is mutated
            genotype = pick_ploids(options.rng, options.ploidy, mut_model.homozygous_freq)

            if location in output_variants_locations:
                # grab the genotype of the existing mutation at that location
                previous_genotype = output_variants[location][1]

                # Try to find a different arrangement. Cap this at 3 tries
                j = 3
                while genotype == previous_genotype or j > 0:
                    options.rng.shuffle(genotype)
                    j -= 1
                if genotype == previous_genotype:
                    continue

            else:
                output_variants[location] = []

            # The logic for this command is find the position to insert with the very fast bisect algorithm,
            # then insert it there. This ensures that the list stays sorted.
            output_variants_locations.insert(bisect.bisect(output_variants_locations, location), location)

            # This will add a second mutation if applicable plus the warning status, or just the mutation.
            output_variants[location].extend([['.', ref, alt, random_mutation_quality_score, random_mutation_filter,
                                             '.', 'GT', '/'.join([str(x) for x in genotype])], genotype])
            # The count will tell us how many we actually added v how many we were trying to add
            n_added += 1

    _LOG.info(f'Finished generating random mutations in {(time.time() - start)/60:.2f} minutes')
    _LOG.info(f'Added {n_added} mutations to {reference.id}')

    _LOG.info(f'Outputting temp vcf for {reference.id} for later use')
    start = time.time()
    filtered_by_target = 0
    filtered_by_discard = 0
    offset = 0

    mutated_references = []
    if options.produce_fasta:
        ref_seq_list = [char for char in reference.seq]
        if options.fasta_per_ploid:
            mutated_references = [ref_seq_list] * options.ploidy
            offset = [0] * options.ploidy
        else:
            mutated_references = ref_seq_list

    # Now that we've generated a list of mutations, let's add them to the temp vcf and fasta (if indicated)
    with open(tmp_vcf_fn, 'a') as tmp:
        for i in range(len(output_variants_locations)):
            variant = output_variants[output_variants_locations[i]]
            # If there is only one variant at this location, this will only execute once. This line skips 2 because in
            # normal vcfs, the dictionary will have one list for te variant, one for the genotype. Will need to adjust
            # for cancer samples
            count = 2
            if options.cancer:
                count = 3

            is_del_list = []
            for j in range(0, len(variant), count):

                # Warn about monomorphic sites but let it go otherwise
                if variant[j][2] == '.':
                    _LOG.warning(f'Found monomorphic reference variant. '
                                 f'These prevent mutations from being added at that location.')

                # Check target and discard bed files
                # Note that if we aren't completely discarding off-target matches, then for the vcf,
                # we don't need to even check if it's in the target region.

                if target_regions and options.discard_offtarget:
                    in_region = False
                    for coords in target_regions:
                        # Check that this location is valid.
                        # Note that we are assuming the coords are half open here.
                        if coords[0] <= output_variants_locations[i] < coords[1]:
                            in_region = True
                            # You can stop looking when you find it
                            break
                    if not in_region:
                        _LOG.debug(f'Variant filtered out by target regions bed: {reference.id}: '
                                   f'{output_variants_locations[i]}')
                        filtered_by_target += 1
                        continue

                # These are firm discards so they all get tossed.
                if discard_regions:
                    in_region = False
                    for coords in discard_regions:
                        # Check if this location is in an excluded zone
                        # Note that we are assuming the coords are half open here
                        if coords[0] <= output_variants_locations[i] < coords[1]:
                            in_region = True
                            break
                    if in_region:
                        _LOG.debug(f'Variant filtered out by discard regions bed: {reference.id}: '
                                   f'{output_variants_locations[i]}')
                        filtered_by_target += 1
                        continue

                if options.produce_fasta:

                    ref = variant[j][1]
                    alt = variant[j][2]
                    alt = alt.split(',')

                    genotype = variant[j + 1]
                    if options.fasta_per_ploid:
                        for k in range(len(genotype)):
                            if genotype[k]:
                                position = output_variants_locations[i] + offset[k]
                                if mutated_references[k][position: position+len(ref)] != ref:
                                    # if our base has been deleted or changed already, we'll skip this one.
                                    continue
                                mutated_references[k] = add_variant_to_fasta(output_variants_locations[i] + offset[k],
                                                                             ref, alt,
                                                                             mutated_references[k])
                                # offset is a running total of the position modification caused by insertions and
                                # deletions. We update it after inserting the variant,
                                # so that the next one is in the correct position.
                                offset[k] += len(alt) - len(ref)
                    else:
                        if 1 in genotype:
                            position = output_variants_locations[i] + offset
                            # If the string of that sequence is different it means our base(s) have been altered.
                            if ''.join(mutated_references[position: position + len(ref)]) != ref:
                                # if our base(s) has been deleted or changed already, we'll skip this one.
                                continue
                            mutated_references = add_variant_to_fasta(output_variants_locations[i] + offset,
                                                                      ref, alt,
                                                                      mutated_references)
                            # offset is a running total of the position modification caused by insertions and deletions
                            # We update it after inserting the variant, so that the next one is in the correct position.
                            offset += len(alt) - len(ref)

                # Add one to get it back into vcf coordinates
                line = f'{reference.id}\t{output_variants_locations[i] + 1}\t{variant[j][0]}\t{variant[j][1]}' \
                       f'\t{variant[j][2]}\t{variant[j][3]}\t{variant[j][4]}\t{variant[j][5]}' \
                       f'\t{variant[j][6]}\t{variant[j][7]}\n'

                if options.cancer:
                    line += f'\t{variant[j][8]}'

                tmp.write(line)

    if filtered_by_target:
        _LOG.info(f'{filtered_by_target} variants excluded because '
                 f'of target regions with discard off-target enabled')

    if filtered_by_discard:
        _LOG.info(f'{filtered_by_discard} variants excluded because '
                  f'of target regions with discard off-target enabled')

    _LOG.info(f'Finished outputting temp vcf in {(time.time() - start)/60:.2f} minutes')

    _LOG.debug(f"Added {n_added} mutations to the reference.")

    # Let's write out the vcf, if asked to
    # TODO move this out
    # if options.produce_vcf:
    #     _LOG.info(f'Writing temp vcf for {reference.id}')
    #
    #     chrom_vcf_file = f"{out_prefix}_{reference.id}.vcf.gz"
    #
    #     # If needed, save a copy of the temp file out as a gzipped vcf.
    #     with open_input(tmp_vcf_fn) as f_in, open_output(chrom_vcf_file) as f_out:
    #         f_out.writelines(f_in)

    return tmp_vcf_fn, mutated_references
