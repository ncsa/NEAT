import gzip
import math
import random
import time
from random import choice
import pathlib
from Bio.Seq import Seq, MutableSeq
import numpy as np
import bisect
import shutil

from source.error_handling import log_mssg, premature_exit
from source.constants_and_defaults import ALLOWED_NUCL
from source.probability import DiscreteDistribution, poisson_list
from source.ploid_functions import which_ploid, pick_ploids
from source.vcf_func import parse_input_vcf
import tempfile
import pandas as pd
import io
import glob
from heapq import merge

import shutil


def close_temp_files(tmp_files):
    for file_handle in tmp_files:
        file_handle.close()


def close_temp_dir(temp_dir):
    temp_dir.cleanup()


def remove_blacklisted(which_ploids, genotype, location, blacklist):
    remove_ploid = []
    if location in blacklist:
        for ploid in which_ploids:
            if ploid in blacklist[location]:
                remove_ploid.append(ploid)
    for ploid in remove_ploid:
        # If there is already a mutation at this location on this ploid, we'll set the genotype
        # for it to zero
        genotype[ploid] = 0
        log_mssg(f'Skipping input variant because a variant '
                 f'already exists at that location {location}, ploid: {ploid}', 'warning')
    return genotype


def parse_mutation_rate_dict(mutation_rate_map, avg_rate, reference):
    """
    This parses the mutation rate dict, in order to fill in the dict so it can by more easily cycled through
    later.

    example mutation_rate_map {'H1N1_HA': [(22, 500, 0.001), (510, 750, 0.003)]}
    example avg_rate: 0.03
    example reference_index: {'H1N1_HA': Seq("AAACA")}
    example non_n_regions: {'H1N1_HA': [(50, 1701)]}

    - intervals that don't overlap
    - intervals with gaps
    - intervals with no gaps
    - mixes of the above
    - empty dictionary ({'chrX': []})

    TODO write tests:
    >>> my_map = [(1, 3, 0.001)]
    >>> my_rate = 0.03
    >>> my_ref = Seq("NAAACAAA")
    >>> x = parse_mutation_rate_dict(my_map, my_rate)
    >>> print(x)
    [(0, 1, 0.03), (1, 3, 0.001), (3, 8, 0.03)]
    """
    # creates the default dict, which is a list of values -1 = not a valid base,
    # any other number is the mutation rate of that base
    ret_list = []
    start = 0
    if mutation_rate_map:
        for region in mutation_rate_map:
            if region[0] > start:
                ret_list.append((start, region[0], avg_rate))
                start = region[1]
                ret_list.append(region)
            elif region[0] == start:
                ret_list.append(region)
                start = region[1]
    else:
        ret_list.append((0, len(reference), avg_rate))
    if ret_list[-1][1] != len(reference):
        ret_list.append((start, len(reference), avg_rate))
    return ret_list


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


def find_random_non_n(mutation_slice, safe_zones, max_attempts):
    for _ in range(max_attempts):
        potential_location = random.randint(0, len(safe_zones) - 1)
        if not safe_zones[potential_location]:
            continue
        else:
            return potential_location + mutation_slice[0]

    return False


def model_trinucs(sequence, models):
    """
    By virtue of the fact that we're trying to speed this up, we'll allow Ns, and
    just set the probs for those at 0
    :param sequence: the sequence to model
    :models: the models for this simulation
    :return: A map of the trinuc probs
    """
    # Set up the model dictionary
    trinuc_models = [0] * len(sequence)
    # We're going to rewrite this to operate on subslices of the reference
    # To do that we need to know the safe zones of just this region.
    # Start at +1 so we get a trinucleotide to start and end one shy for the same reason
    for i in range(1, len(sequence) - 1):
        trinuc = sequence[i - 1:i + 2].seq
        # Let's double check to make sure we didn't pick up a stray N
        if any([j for j in trinuc if j not in ALLOWED_NUCL]):
            trinuc_models[i] = 0.0
        trinuc_models[i] = models.mutation_model['trinuc_bias'][trinuc]
    trinuc_models = np.array(trinuc_models)
    # What if there are valid bases, but no valid trinculeotides? Skip this one.
    if not any(trinuc_models):
        return None
    # else, make a discrete distribution
    return DiscreteDistribution(range(len(sequence)), trinuc_models)


def split_sequence(my_string):
    return [char for char in my_string]


def check_if_deleted(all_dels, location):

    for deletion in all_dels:
        # deletion[0] = location of deletion
        # deletion[1] = length of the deletion
        if deletion[0] < location < deletion[0] + deletion[1]:
            return True

    return False


def generate_variants(reference, chrom, tmp_vcf_fn, target_regions, discard_regions,
                      mutation_rate_regions, output_variants, models, options,
                      out_prefix):
    """
    This function will take all the setup we did in the main part and actually do NEAT stuff.

    TODO: need to add cancer logic to this section
    """

    # Step 1: Create a VCF of mutations
    log_mssg(f'Generating chromosome map', 'info')
    mutation_map = parse_mutation_rate_dict(mutation_rate_regions, models.mutation_model['avg_mut_rate'], reference)

    mutation_regions_model = DiscreteDistribution([(x[0], x[1]) for x in mutation_map],
                                                  np.array([x[2] for x in mutation_map]))

    # Trying to use a random window to keep memory under control. May need to adjust this number.
    max_window_size = 1000
    n_added = 0
    random_mutation_quality_score = max(models.sequencing_error_model.quality_scores)
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

    log_mssg(f'Adding random mutations for {chrom}', 'info')
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
    factors.append(total * models.mutation_model['avg_mut_rate'])
    overall_mutation_average = sum(factors)/len(reference)
    average_number_of_mutations = int(len(reference) * overall_mutation_average)
    max_mutations = int(len(reference) * 0.3)
    total_mutations_model = poisson_list(max_mutations, average_number_of_mutations)

    # This number will serve as a counter for our loops. Let's add a minimum of 1 mutation no matter what.
    min_mutations = 1
    if options.min_mutations:
        min_mutations = options.min_mutations

    how_many_mutations = total_mutations_model.sample() + min_mutations

    log_mssg(f'Planning to add {how_many_mutations} mutations. The final number may be less.', 'info')

    # We may need to skip locations if they were deleted
    all_dels = []

    while how_many_mutations > 0:
        # Pick a region based on the mutation rates
        # (default is one rate for the whole chromosome, so this will be trivial in that case
        mut_region = mutation_regions_model.sample()
        # Pick a random starting place. Randint is inclusive of endpoints, so we subtract 1
        window_start = random.randint(mut_region[0], mut_region[1] - 1)
        found = False
        if reference[window_start] not in ALLOWED_NUCL:
            # pick a random location to the right
            plus = random.randint(window_start + 1, mut_region[1] - 1)
            if reference[plus] in ALLOWED_NUCL:
                found = True
                window_start = plus
            else:
                # If that didn't work pick a random location to the left
                if window_start - 1 > mut_region[0]:
                    minus = random.randint(mut_region[0], window_start - 1)
                    if reference[minus] in ALLOWED_NUCL:
                        found = True
                        window_start = minus
        else:
            found = True

        # If we couldn't find a spot, try again
        if not found:
            log_mssg(f"Couldn't find a spot for this one", 'debug')
            continue

        # at this point the location is found. Now we need a second endpoint. Grab 1000 bases or to the end.
        end_point = random.randint(window_start, min(mut_region[1], window_start + max_window_size) - 1)
        if reference[end_point] not in ALLOWED_NUCL:
            # Didn't find it to the right, look to the left
            end_point = random.randint(max(window_start - max_window_size, mut_region[0]), window_start)
            if reference[end_point] not in ALLOWED_NUCL:
                # No suitable end_point, so we try again
                log_mssg(f"No suitable end_point", 'debug')
                continue

        # Sorting assures that wherever we found the end point, the coordinates will be in the correct order for slicing
        mutation_slice = sorted([window_start, end_point])
        slice_distance = mutation_slice[1] - mutation_slice[0]
        # How many variants to add in this slice (at least one, or we'll never get to the finish line)
        variants_to_add_in_slice = max(int((slice_distance/len(reference)) * how_many_mutations), 1)
        # log_mssg(f"Planning to add {variants_to_add_in_slice} variants to slice", 'debug')

        subsequence = reference[mutation_slice[0]: mutation_slice[1]]
        # In the case where the end points are equal, just skip
        if len(subsequence) == 0:
            log_mssg(f"Subsequence had a zero length", 'debug')
            continue

        non_n_regions = map_non_n_regions(subsequence)

        # If the sequence has too many N's, we'll skip it
        if not non_n_regions:
            log_mssg(f"In N region", 'debug')
            continue

        # Begin random mutations for this slice
        while variants_to_add_in_slice > 0:
            # We decrement now because we don't want to get stuck in a never ending loop
            variants_to_add_in_slice -= 1
            how_many_mutations -= 1

            # Now figure out the type of random mutation to insert
            is_indel = random.random() <= models.mutation_model['indel_freq']
            # Case 1: indel
            if is_indel:
                # First pick a location. This function ensures that we do not pick an N as our starting place
                # Note that find_random_non_n helpfully adds the slice start to the location.
                location = find_random_non_n(mutation_slice, non_n_regions, 5)

                if check_if_deleted(all_dels, location):
                    continue  # No increments, no attempts, just try again.

                is_insertion = random.random() <= models.mutation_model['indel_insert_percentage']
                if is_insertion:
                    length = models.mutation_model['insert_length_model'].sample()
                    # Try to find a location to insert. Give it ten tries, then give up.
                    ref = reference[location]
                    insertion = ''.join(random.choices(ALLOWED_NUCL, k=length))
                    alt = ref + insertion
                else:
                    length = models.mutation_model['deletion_length_model'].sample()
                    # Plus one so we make sure to grab the first base too
                    ref = reference[location: location+length+1].seq
                    alt = reference[location]
                    all_dels.append((location, length))

            # Case 2: SNP
            else:
                trinuc_probs = model_trinucs(subsequence, models)
                # If we have some edge case where there was no actual valid trinucleotides, we'll skip this
                if not trinuc_probs:
                    log_mssg(f"Could not build trinuc probs", 'debug')
                    continue
                # So there is at least valid trinuc in this subsequence, so now we sample for a location
                # It's a relative location, so we add the start point of the subsequence to that.
                location = trinuc_probs.sample() + mutation_slice[0]

                if check_if_deleted(all_dels, location):
                    continue  # No increments, no attempts, just try again.

                ref = reference[location]
                trinuc = reference[location - 1: location + 2]
                transition_values = ALLOWED_NUCL
                transition_probs = models.mutation_model['trinuc_trans_prob'][trinuc.seq].values()
                # We want this to be an actual variant, so we'll try a few times
                alt = random.choices(transition_values, transition_probs)[0]
                # Max 10 tries
                j = 10
                while alt == ref and j > 0:
                    alt = random.choices(transition_values, transition_probs)[0]
                    j -= 1

                # We tried and failed, so we give up.
                if alt == ref:
                    log_mssg(f"Transition failed to find an alt", 'debug')
                    continue

            # Now we have a location, an alt and a ref, so we figure out which ploid is mutated
            genotype = pick_ploids(options.ploidy, models.mutation_model['homozygous_freq'])

            if location in output_variants_locations:
                # grab the genotype of the existing mutation at that location
                previous_genotype = output_variants[location][1]

                # Try to find a different arrangement. Cap this at 3 tries
                j = 3
                while genotype == previous_genotype or j > 0:
                    random.shuffle(genotype)
                    j -= 1
                if genotype == previous_genotype:
                    log_mssg(f"Skipping random mutation because a variant already exists there", 'debug')
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

    log_mssg(f'Finished generating random mutations in {(time.time() - start)/60:.2f} minutes', 'info')
    log_mssg(f'Added {n_added} mutations to {chrom}', 'info')

    log_mssg(f'Outputting temp vcf for {chrom} for later use', 'info')
    start = time.time()
    filtered_by_target = 0
    filtered_by_discard = 0
    offset = 0

    mutated_references = []
    if options.produce_fasta:
        ref_seq_list = split_sequence(str(reference.seq))
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
                    log_mssg(f'Found monomorphic reference variant. '
                             f'These prevent mutations from being added at that location.', 'info')

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
                        log_mssg(f'Variant filtered out by target regions bed: {chrom}: '
                                 f'{output_variants_locations[i]}', 'debug')
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
                        log_mssg(f'Variant filtered out by discard regions bed: {chrom}: '
                                 f'{output_variants_locations[i]}', 'debug')
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
                                # offset is a running total of the position modification caused by insertions and deletions
                                # We update it after inserting the variant, so that the next one is in the correct position.
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
                line = f'{chrom}\t{output_variants_locations[i] + 1}\t{variant[j][0]}\t{variant[j][1]}' \
                       f'\t{variant[j][2]}\t{variant[j][3]}\t{variant[j][4]}\t{variant[j][5]}' \
                       f'\t{variant[j][6]}\t{variant[j][7]}\n'

                if options.cancer:
                    line += f'\t{variant[j][8]}'

                tmp.write(line)

    if filtered_by_target:
        log_mssg(f'{filtered_by_target} variants excluded because '
                 f'of target regions with discard off-target enabled', 'info')

    if filtered_by_discard:
        log_mssg(f'{filtered_by_discard} variants excluded because '
                 f'of target regions with discard off-target enabled', 'info')

    log_mssg(f'Finished outputting temp vcf in {(time.time() - start)/60:.2f} minutes', 'info')

    log_mssg(f"Added {n_added} mutations to the reference.", 'debug')

    # Let's write out the vcf, if asked to
    if options.produce_vcf:
        log_mssg(f'Creating vcf for {chrom}', 'info')

        chrom_vcf_file = f"{out_prefix}_{chrom}.vcf.gz"

        # If needed, save a copy of the temp file out as a gzipped vcf.
        with open(tmp_vcf_fn, 'rb') as f_in, gzip.open(chrom_vcf_file, 'wb') as f_out:
            f_out.writelines(f_in)

        return chrom_vcf_file, tmp_vcf_fn, mutated_references

    return None, tmp_vcf_fn, mutated_references


def add_variant_to_fasta(position: int, ref: str, alt: str, reference_to_alter: list):
    ref_len = len(ref)
    # Apply mutation
    reference_to_alter[position:position+ref_len] = alt
    return reference_to_alter


def generate_reads(reference_chrom, models, input_vcf, temporary_directory, options, chrom):
    # apply errors and mutations and generate reads
    chrom_fastq_r1 = temporary_directory / f'{chrom}_tmp_r1.fq'
    chrom_fastq_r2 = None
    if options.paired_ended:
        chrom_fastq_r2 = temporary_directory / f'{chrom}_tmp_r2.fq'

    # # I think I need all the non-n regions for this analysis
    # non_n_regions = map_non_n_regions(reference_chrom)
    # total_bp_spon = sum(non_n_regions)
    #
    # current_progress = 0
    # current_percent = 0
    # have_printed_100 = False
    #
    # all_variants_out = {}
    # sequences = None
    # window_target_scale = 100
    #
    # target_size = window_target_scale * options.read_len
    # overlap = options.read_len
    # overlap_min_window_size = options.read_len + 10
    #
    # if options.paired_ended:
    #     target_size = window_target_scale * options.fragment_mean
    #     overlap = options.fragment_mean
    #     overlap_min_window_size = max(models.fraglen_model.values) + 10

    scale_fator = options.fragment_mean if options.paired_ended else options.read_len
    window_scale = 100 * scale_fator

    vars_to_insert = {}
    with gzip.open(input_vcf, 'r') as input_variants:
        for line in input_variants:
            if line.startswith("@") or line.startswith("#"):
                continue
            line_split = line.strip().split('\t')
            # Since these vars are on the same chromosome, we can index by position
            # We don't need the ID for this so let's skip it
            # The remaining fields will be:
            #   - vars_to_insert[position][0] = REF
            #   - vars_to_insert[position][1] = ALT
            #   - vars_to_insert[position][2] = QUAL
            #   - vars_to_insert[position][3] = FILTER
            #   - vars_to_insert[position][4] = INFO
            #   - vars_to_insert[position][5] = FORMAT
            #   - vars_to_insert[position][6] = SAMPLE_1
            #   - vars_to_insert[position][7] = SAMPLE_2 (optional)
            vars_to_insert[line_split[1]] = line_split[3:]

    log_mssg(f'Sampling reads...', 'info')
    start = time.time()
    """
    I'm going to try to use a percentage N for this. Basically, if there are at least 40% 
    Real reads, then we'll use it, otherwise we won't use it. I have no idea how illumina actually works
    in this regard.
    """
    i = 0
    while i < len(reference_chrom):
        # window = reference_chrom[i: i + target_size]
        # # calculate the percentage of valid nucleotides in the window, if it's less than 40%, skip this window
        # if sum(window)/len(window) < 0.4:
        #     i += 1
        #     continue
        # sample_reads(window)
        if reference_chrom[i] == "N":
            i += 1
            continue

        # Found a non-N
        start = i - random.randint(0, 25)
        end = start + window_scale

        i += 1

    log_mssg(f"Finished sampling reads in {time.time() - start} seconds", 'info')
    return chrom_fastq_r1, chrom_fastq_r2


def generate_bam(reference_chrom, options, temporary_vcf, temporary_dir, chrom):
    # Step 4 (optional): Create a golden bam with the reads aligned to the original reference
    pass

