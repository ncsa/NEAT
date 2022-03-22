import gzip
import math
import random
import time
from random import choice
import pathlib
from Bio.Seq import Seq
import numpy as np
import bisect

from source.error_handling import log_mssg, premature_exit
from source.constants_and_defaults import ALLOWED_NUCL
from source.probability import DiscreteDistribution, poisson_list
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


def pick_ploids(ploidy, homozygous_frequency, number_alts=1) -> list:
    """
    Applies alts to random ploids. Picks at least one, maybe more.
    :param ploidy: how many copies of each chromosome this organism has
    :param number_alts: If there is more than one alt, this will assign a random alt to the ploids
    :return: a list of strings representing the genotype of each ploid.
    """
    # number of ploids to make this mutation on (always at least 1)
    if random.random() < homozygous_frequency:
        if ploidy <= 2:
            # If it's homozygous. it's on both ploids
            how_many = ploidy
        else:
            # if it's polyploid, we'll just pick some.
            # TODO may need to improve the modeling for polyploid
            how_many = 1
            for i in range(ploidy):
                # Not totally sure how to model this, so I'm counting each
                # ploid as a separate homozygous event. That doesn't exactly make
                # sense though so we'll improve this later.
                if random.random() < homozygous_frequency:
                    how_many += 1
                else:
                    break
    else:
        how_many = 1

    # wp is just the temporary genotype list
    wp = [0] * ploidy
    while how_many > 0:
        x = random.choice(range(ploidy))
        # pick a random alt. in VCF terminology, 0 = REF, 1 = ALT1, 2 = ALT2, etc
        wp[x] = random.choice(range(1, number_alts + 1))
        how_many -= 1

    return [str(x) for x in wp]


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


def find_random_non_n(slice, safe_zones, max_attempts):
    for _ in range(max_attempts):
        potential_location = random.randint(slice[0], slice[1] - 1)
        if not safe_zones[potential_location]:
            continue
        else:
            return potential_location

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


def execute_neat(reference, chrom, out_prefix_name, target_regions, discard_regions,
                 mutation_rate_regions, input_variants, models, options,
                 out_prefix):
    """
    This function will take all the setup we did in the main part and actually do NEAT stuff.
    """
    # Might be able to pre-populate this
    final_files = []

    # Setting up temp files to write to
    tmp_fasta_fn = None
    tmp_fastq1_fn = None
    tmp_fastq2_fn = None
    tmp_sam_fn = None
    tmp_vcf_fn = None

    # Since we're only running single threaded for now:
    threadidx = 1
    # rough count of the number of batches we want to make for this
    num_batches = len(reference)//1000000 + 1

    temporary_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = pathlib.Path(temporary_dir.name).resolve()

    # We need the temp vcf file no matter what
    tmp_vcf_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{chrom}_{threadidx}.vcf"
    log_mssg(f'tmp_vcf_fn = {tmp_vcf_fn}', 'debug')

    if options.produce_bam:
        tmp_sam_fn = tmp_dir_path / f"{out_prefix_name}_tmp_records_{threadidx}.tsam"
        log_mssg(f'tmp_sam_fn = {tmp_sam_fn}', 'debug')
    if options.produce_fasta:
        tmp_fasta_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}.fasta"
        log_mssg(f'tmp_fasta_fn = {tmp_fasta_fn}', 'debug')
    if options.produce_fastq:
        tmp_fastq1_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}_read1.fq"
        log_mssg(f'tmp_fastq1_fn = {tmp_fastq1_fn}', 'debug')
        if options.paired_ended:
            tmp_fastq2_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}_read2.fq"
            log_mssg(f'tmp_fastq2_fn  = {tmp_fastq2_fn}', 'debug')

    if threadidx == 1:
        # init_progress_info()
        pass

    # Step 1: Create a VCF of variants (mutation and sequencing error variants)
    # We'll create a temp file first then output it if the user requested the file
    with open(tmp_vcf_fn, 'w') as tmp_vcf_file:
        tmp_vcf_file.write(f'@NEAT temporary file, for generating the list of mutations.\n')

    log_mssg(f'Generating chromosome map', 'info')
    mutation_map = parse_mutation_rate_dict(mutation_rate_regions, models.mutation_model['avg_mut_rate'], reference)

    mutation_regions_model = DiscreteDistribution([(x[0], x[1]) for x in mutation_map],
                                                  np.array([x[2] for x in mutation_map]))

    start = time.time()

    # Trying te prevent inserting a mutation in more than one place
    blacklist = {}
    max_window_size = 1000
    n_added = 0
    random_mutation_quality_score = max(models.sequencing_error_model.quality_scores)
    random_mutation_filter = "PASS"
    """
    Inserting any input mutations
    """
    if input_variants:
        """
        For reference, the columns in input_variants, and their indices:
        key: (CHROM, POS)
        values[0]: ID [0], REF [1], ALT [2], QUAL [3], FILTER [4], INFO [5], 
                FORMAT [6, optional], SAMPLE1 [7, optional], SAMPLE2 [8, optional]
        values[1]: genotype (if present in original vcf), None if not present
        values[2]: genotype tumor (if present in original vcf), None if not present and a cancer sample
        """
        # this should give is just the positions of the inserts.
        input_variants_locations = input_variants.keys()

    log_mssg(f'Adding mutations for {chrom}', 'info')
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
    how_many_mutations = total_mutations_model.sample()

    log_mssg(f'Mutating chrom: {chrom}...', 'info')

    # Initialize to include all input variants
    all_mutations = input_variants
    while how_many_mutations > 0:
        # Pick a region based on the mutation rates
        # (default is one rate for the whole chromosome, so this will be trivial in that case
        mut_region = mutation_regions_model.sample()
        # Pick a random starting place. Randint is inclusive of endpoints, so we subtract 1
        window_start = random.randint(mut_region[0], mut_region[1] - 1)
        found = False
        if reference[window_start] not in ALLOWED_NUCL:
            # pick a random location to the right
            plus = random.randint(window_start, mut_region[1] - 1)
            if reference[plus] in ALLOWED_NUCL:
                found = True
                window_start = plus
            else:
                # If that didn't work pick a random location to the left
                minus = random.randint(mut_region[0], window_start)
                if reference[minus] in ALLOWED_NUCL:
                    found = True
                    window_start = minus
        else:
            found = True

        # If we couldn't find a spot, try again
        if not found:
            continue

        # at this point the location is found. Now we need a second endpoint. Grab a 1000 bases or to the end.
        end_point = random.randint(window_start, min(mut_region[1], window_start + max_window_size) - 1)
        if reference[end_point] not in ALLOWED_NUCL:
            # Didn't find it to the right, look to the left
            end_point = random.randint(max(window_start - max_window_size, mut_region[0]), window_start - 1)
            if reference[end_point] not in ALLOWED_NUCL:
                # No suitable end_point, so we try again
                continue

        # Sorting assures that wherever we found the end point, the coordinates will be in the correct order for slicing
        mutation_slice = sorted([window_start, end_point])
        slice_distance = mutation_slice[1] - mutation_slice[0]
        # How many variants to add in this slice (at least one, or we'll never get to the finish line)
        variants_to_add_in_slice = max((slice_distance/len(reference)) * how_many_mutations, 1)

        subsequence = reference[mutation_slice[0]: mutation_slice[1]]
        # In the case where the end points are equal, just skip
        if len(subsequence) == 0:
            continue

        non_n_regions = map_non_n_regions(subsequence)

        # If the sequence has too many N's, we'll skip it
        if not non_n_regions:
            continue

        # First find any input variants that go in this window, initialize with input variants
        variants_to_add = []
        if input_variants:
            variants_to_add.append(sorted([x for x in input_variants_locations
                                   if mutation_slice[0] <= x < mutation_slice[1]]))

        # Begin random mutations
        for _ in range(variants_to_add_in_slice):
            # Now figure out the type of random mutation to insert
            is_indel = random.random() <= models.mutation_model['indel_freq']
            is_insertion = False
            length = 0
            # Case 1: indel
            if is_indel:
                # First pick a location. This function ensures that we do not pick an N as our starting place
                location = find_random_non_n(mutation_slice, non_n_regions, 10)
                if not location:
                    continue

                is_insertion = random.random() <= models.mutation_model['indel_insert_percentage']
                if is_insertion:
                    length = models.mutation_model['insert_length_model'].sample()
                    # Try to find a location to insert. Give it ten tries, then give up.
                    ref = reference[location]
                    insertion = ''.join(random.choices(ALLOWED_NUCL, k=length))
                    alt = ref + insertion
                else:
                    length = models.mutation_model['deletion_length_model'].sample()
                    # Now we know what it is and how long, we just need to find a spot for it.
                    ref = ""
                    for i in range(length):
                        # At this point if we find an N, replace with random base and schedule for deletion
                        if not non_n_regions[location + i]:
                            ref += random.choice(ALLOWED_NUCL)
                        else:
                            ref += reference[location + i]
                    alt = reference[location]

            # Case 2: SNP
            else:
                trinuc_probs = model_trinucs(subsequence, models)
                # If we have some edge case where there was no actual valid trinucleotides, we'll skip this
                if not trinuc_probs:
                    continue
                # So there is at least valid trinuc in this subsequence, so now we sample for a location
                # It's a relative location, so we add the start point of the subsequence to that.
                location = trinuc_probs.sample() + mutation_slice[0]
                ref = reference[location]
                trinuc = reference[location - 1: location + 2]
                transition_values = ALLOWED_NUCL
                transition_probs = models.mutation_model['trinuc_trans_prob'][trinuc.seq].values()
                # We want this to be an actual variant, so we'll try a few times
                for _ in range(10):
                    alt = random.choices(transition_values, transition_probs)[0]
                    if alt != ref:
                        break
                # We tried and failed, so we give up.
                if alt == ref:
                    continue
            warning = False
            if location in variants_to_add:
                warning = True

            # Now we put it all together for this window
            random_ploid = False

            # TODO we'll have to add extra code for cancer, later
            # If there are no input variants, this list is only 1 element long.
            # The logic for this command is find the position to insert with the very fast bisect algorithm,
            # then insert it there. This ensures that the list stays sorted.
            variants_to_add.insert(bisect.bisect(variants_to_add, location), location)
            # If there was already a variant at that location, need to pick a ploid that is clear
            if warning:
                # grab the genotype of the existing mutation at that location
                input_genotype = input_variants[location][1]
                if input_genotype:
                    # Check which ploids have mutations available
                    which_ploid = [input_genotype.index(x) for x in input_genotype if x == 0]
                    if which_ploid:
                        ploid = random.choice(which_ploid)
                        genotype = [0] * options.ploidy
                        genotype[ploid] = 1
                        genotype = [str(x) for x in genotype]
                    else:
                        random_ploid = True
            else:
                random_ploid = True

            # Otherwise, determine genotype randomly
            if random_ploid:
                genotype = pick_ploids(options.ploidy, models.mutation_model['homozygous_freq'])

            if location not in all_mutations:
                all_mutations[location] = []
            # This will add a second mutation if applicable plus the warning status, or just the mutation.
            all_mutations[location].extend([['.', ref, alt, random_mutation_quality_score, random_mutation_filter,
                                             '.', 'GT', ''.join(genotype)], genotype])

        # Now that we've generated a list of mutations, let's add them to the temp vcf
        # TODO add cancer logic
        for position in variants_to_add:
            variant = all_mutations[position]

            # If there is only one variant at this location, this will only execute once
            for i in range(0, len(variant), 2):
                ref_sequence = reference[position:position + len(variant[i][1])].seq
                # Some checks to make sure we are inserting into a good spot

                # Check if we're somewhere we shouldn't be
                allowed = reference[position] in ALLOWED_NUCL

                # Make sure the ref field is correct
                if ref_sequence != variant[i][1]:
                    log_mssg(f"Skipping variant where reference does not match "
                             f"input vcf at {chrom}: {position} - {variant[0][1]} "
                             f"(ref = {ref_sequence}", 'warning')
                    continue

                # Check target and discard bed files
                # Note that if we aren't completely discarding off-target matches, then for the vcf,
                # we don't need to even check if it's in the target region.
                if target_regions and options.discard_offtarget and allowed:
                    in_region = False
                    for coords in target_regions:
                        # Check that this location is valid.
                        # Note that we are assuming the coords are half open here.
                        if coords[0] <= position < coords[1]:
                            in_region = True
                            # You can stop looking when you find it
                            break
                    if not in_region:
                        allowed = False

                # These are firm discards so they all get tossed.
                if discard_regions and allowed:
                    in_region = False
                    for coords in discard_regions:
                        # Check if this location is in an excluded zone
                        # Note that we are assuming the coords are half open here
                        if coords[0] <= position < coords[1]:
                            in_region = True
                            break
                    if in_region:
                        allowed = False

                # At this point, if it's not allowed, just skip
                if not allowed:
                    continue

                if variant[i][2] == '.':
                    log_mssg(f'Found monomorphic reference variant. '
                             f'These prevent mutations from being added at that location.', 'info')
                else:
                    number_of_alts = len(variant[i][2].split(','))
                    # If genotype data was not present in the vcf
                    # We already have genotype data for input vcfs, so this will handle that case too
                    if not variant[i+1]:
                        genotype = pick_ploids(options.ploidy, number_of_alts)
                    else:
                        genotype = variant[i+1]
                    if options.cancer:
                        if not variant[2]:
                            genotype_tumor = pick_ploids(options.ploidy, number_of_alts)
                        else:
                            genotype_tumor = variant[i+2]

                # determine which ploid was mutated, for the blacklist. Find all the mutated ploids (genotype != 0)
                which_ploids = [genotype.index(x) for x in genotype if x != 0]
                if options.cancer:
                    tumor_ploids = [genotype_tumor.index(x) for x in genotype_tumor if x != 0]

                genotype = remove_blacklisted(which_ploids, genotype, position, blacklist)

                # if that leaves us with all zeros, skip this one
                if not any(genotype):
                    continue
                else:
                    # Which ploids are still being mutated?
                    which_ploids = [genotype.index(x) for x in genotype if x != 0]
                # That is, if we have format and sample columns, we need to add in the GT
                if len(variant[i]) > 6 and not options.cancer:
                    if "GT" not in variant[i][6]:
                        format_field = f"GT:{variant[i][6]}"
                        gt_index = 0
                    else:
                        format_field = variant[i][6]
                        gt_index = variant[i][6].split(':').index("GT")
                    # Slice up the field and add in the genotype, then join it back together
                    sample_field = variant[i][7].split(':')
                    gt_field = sample_field[:gt_index] + ["/".join(genotype)] + sample_field[gt_index:]
                    gt_field = ":".join(gt_field)
                else:
                    format_field = variant[0][6]
                    gt_field = variant[0][7]

                tail = f'{format_field}\t{gt_field}'

                line = f'{chrom}\t{position}\t{variant[i][0]}\t{variant[i][1]}\t' \
                       f'{variant[i][2]}\t{variant[i][3]}\t{variant[i][4]}\t{variant[i][5]}\t' \
                       f'{tail}\n'
                how_many_mutations -= 1

                with open(tmp_vcf_fn, 'a') as tmp:
                    tmp.write(line)

                # Initialize blacklist for that position
                if position not in blacklist:
                    blacklist[position] = []
                    for index in which_ploids:
                        blacklist[position].append(index)


    log_mssg(f"Finished mutating {chrom}. Time: {time.time() - start}", 'debug')
    log_mssg(f"Added {n_added} mutations to the reference.", 'debug')

    # Let's write out the vcf, if asked to
    if options.produce_vcf:
        log_mssg(f'Sorting and writing temp vcf for {chrom}', 'info')

        path = f"{tmp_dir_path}/chunk_*.vcf"
        chunksize = 1_000_000
        fid = 1
        lines = []
        number_of_files = math.ceil(len(all_mutations) / chunksize)

        with open(str(tmp_vcf_fn), 'r') as f:
            f_out = open(f'{tmp_dir_path}/chunk_{fid}.vcf', 'w')
            for line_num, line in enumerate(f, 1):
                if not line.startswith('@'):
                    lines.append(line)
                if not line_num % chunksize:
                    log_mssg(f"Sorting file chunk {fid} of {number_of_files}", 'debug')
                    lines = sorted(lines, key=lambda k: (k.split()[0], int(k.split()[1])))
                    f_out.writelines(lines)
                    f_out.close()
                    lines = []
                    fid += 1
                    f_out = open(f'{tmp_dir_path}/chunk_{fid}.vcf', 'w')
            # If the above loop has a number that is not a multiple of chunksize, then it will
            # have a remainder. In that case, we do the final step one last time:
            if lines:
                log_mssg(f"Sorting file chunk {fid} {number_of_files}", 'debug')
                lines = sorted(lines, key=lambda k: (k.split()[0], int(k.split()[1])))
                f_out.writelines(lines)
                f_out.close()

        chunks = []
        for filename in glob.glob(path):
            chunks += [open(filename, 'r')]

        log_mssg(f'Sorting chunks into vcf file for {chrom}', 'debug')
        file = f"{out_prefix}_{chrom}.vcf.gz"
        with gzip.open(file, 'wt') as f_out:
            f_out.writelines(merge(*chunks, key=lambda k: (k.split()[0], int(k.split()[1]))))

        for item in chunks:
            item.close()

    temporary_dir.cleanup()

    # Step 2 (optional): create a fasta with those variants inserted
    # Step 3 (optional): create a fastq from the mutated fasta
    # Step 4 (optional): Create a golden bam with the reads aligned to the original reference


    # TODO add multithreading
    # with WorkerPool(n_jobs=options.threads, shared_objects=jobs) as pool:
    #     pool.map(run_neat, jobs, progress_bar=True)

    return file
