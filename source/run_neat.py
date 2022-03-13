import gzip
import random
import time
from random import choice
import pathlib
import numpy as np
from mpire import WorkerPool

from source.error_handling import print_and_log, premature_exit
from source.constants_and_defaults import ALLOWED_NUCL
from source.probability import DiscreteDistribution
from source.vcf_func import parse_vcf
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


def pick_ploids(ploidy, homozygous_freq, number_alts=1) -> list:
    """
    Applies alts to random ploids. Picks at least one, maybe more.
    :param ploidy: how many copies of each chromosome this organism has
    :param number_alts: If there is more than one alt, this will assign a random alt to the ploids
    :return: a list of strings representing the genotype of each ploid.
    """
    # number of ploids to make this mutation on (always at least 1)
    if random.random() < homozygous_freq:
        # If it's homozygous. it's on all the ploids
        # TODO fact check this with Christina
        how_many = ploidy
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

    temporary_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = pathlib.Path(temporary_dir.name).resolve()

    # We need the temp vcf file no matter what
    tmp_vcf_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{chrom}_{threadidx}.vcf"
    if options.debug:
        print_and_log(f'tmp_vcf_fn = {tmp_vcf_fn}', 'debug')

    if options.produce_bam:
        tmp_sam_fn = tmp_dir_path / f"{out_prefix_name}_tmp_records_{threadidx}.tsam"
        if options.debug:
            print_and_log(f'tmp_sam_fn = {tmp_sam_fn}', 'debug')
    if options.produce_fasta:
        tmp_fasta_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}.fasta"
        if options.debug:
            print_and_log(f'tmp_fasta_fn = {tmp_fasta_fn}', 'debug')
    if options.produce_fastq:
        tmp_fastq1_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}_read1.fq"
        if options.debug:
            print_and_log(f'tmp_fastq1_fn = {tmp_fastq1_fn}', 'debug')
        if options.paired_ended:
            tmp_fastq2_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}_read2.fq"
            if options.debug:
                print_and_log(f'tmp_fastq2_fn  = {tmp_fastq2_fn}', 'debug')

    if threadidx == 1:
        # init_progress_info()
        pass


    # Step 1: Create a VCF of variants (mutation and sequencing error variants)
    # We'll create a temp file first then output it if the user requested the file
    with open(tmp_vcf_fn, 'w') as tmp_vcf_file:
        tmp_vcf_file.write(f'@NEAT temporary file, for generating the list of mutations.\n')

    contig_sequence = reference.seq

    start = time.time()
    """
    Inserting any input mutations
    """
    for variant in input_variants.iterrows():
        ref_sequence = contig_sequence[variant.POS:len(variant.REF)]
        if ref_sequence != variant.REF:
            print_and_log(f"Skipping variant where reference does not match "
                          f"input vcf at {chrom}: {variant.POS}", 'warning')
            continue
        number_of_alts = len(variant.ALT.split(','))
        genotype = pick_ploids(options.ploidy, number_of_alts)

        line = f'{chrom}\t{variant.POS}\t{variant.ID}\t{variant.REF}\t' \
               f'{variant.ALT}\t{variant.QUAL}\t' \
               f'PASS\t{variant.INFO}\t' \
               f'GT\t{"/".join(genotype)}\n'

        with open(tmp_vcf_fn, 'a') as tmp:
            tmp.write(line)

    if options.debug:
        print_and_log(f'Completed inserting the input mutations for {chrom} in {time.time() - start}', 'debug')

    print_and_log(f'Mutating {chrom}', 'info')
    start = time.time()

    """
    Decide how many and where any random mutations will happen in this contig
    """
    # Create a dictionary of regions. With no mutation_rate_regions, then there is only one region for the contig
    # When we add multithreading, I think we can use these intervals for parallel processing.
    # Establish a uniform mutation rate by default
    mutation_regions = [(0, len(contig_sequence), models.mutation_model['avg_mut_rate'])]
    overall_mutation_rate = models.mutation_model['avg_mut_rate']
    if mutation_rate_regions:
        mutation_regions = mutation_rate_regions
        # This uses an average across the contig to calculate how many total variants to add. At least 1.
        overall_mutation_rate = sum([x[2] for x in mutation_regions])/len(mutation_regions)

    # Figure out the weight of each region. If there is only one region, then this will be trivial
    weighted_mutation_regions = {x: (x[1] - x[0]) * x[2] for x in mutation_regions}
    # This will sort this dict by weight
    weighted_mutation_regions = dict(sorted(weighted_mutation_regions.items(), key=lambda x: x[1]))

    # This model will allow us to figure out what regions to place mutations in
    mutation_regions_model = DiscreteDistribution(list(weighted_mutation_regions.keys()),
                                                  weighted_mutation_regions.values())

    mutations_to_add = int(len(contig_sequence) * overall_mutation_rate) + 1
    if options.debug:
        print_and_log(f'Planning to add {mutations_to_add} mutations to {chrom}', 'debug')
    
    mutation_data = {x: [] for x in range(mutations_to_add)}

    print_and_log(f'Generating mutation positions.', 'info')
    for variant in range(mutations_to_add):
        genotype = pick_ploids(options.ploidy, models.mutation_model['homozygous_freq'])
        region = mutation_regions_model.sample()
        # for now our options are indel or snp. Later we can add more variants.
        is_indel = random.random() <= models.mutation_model['indel_freq']
        is_insertion = False
        length = 0
        if is_indel:
            is_insertion = random.random() <= models.mutation_model['indel_insert_percentage']
            if is_insertion:
                length = models.mutation_model['insert_length_model'].sample()
            else:
                length = models.mutation_model['deletion_length_model'].sample()
        # Find somewhere to put this damn thing
        if is_indel:
            # Indels can go anywhere in the region
            potential_location = random.randint(range(region[0], region[1]))
        else:
            # Use the trinuc bias to find a spot for this SNP

            # where to look
            subsequence = contig_sequence[region[0]:region[1]]
            # try 100 times to find an appropriate trinuc
            found = False
            seek = -1
            max_attempts = 100
            while max_attempts > 0:
                max_attempts -= 1
                trinuc_to_mutate = models.mutation_model['trinuc_mut_prob'].sample()
                # This index will be relative to the start of the region,
                # so we add on the region start to keep the coords relative to the reference.
                seek = subsequence.find(trinuc_to_mutate) + region[0]
                if seek == -1:
                    continue
                else:
                    # At least on occurrence of this trinucleotide was found in this set, so
                    # we'll have to pick one at random
                    total = subsequence.count_overlap(trinuc_to_mutate)
                    pick_one = 0
                    if total > 1:
                        pick_one = random.randint(1, total)
                    for i in range(pick_one):
                        # seek is set to the first instance. If we picked 0, then
                        # that will be our choice and this loop will skippde
                        # otherwise, we check subsequences until we find the correct trinuc
                        # to use. We add one to the index to skip the current trinuc
                        subsubsequence = subsequence[seek + 1:]
                        # Add one to match the subsequence check in the previous line
                        seek += subsubsequence.find(trinuc_to_mutate) + 1
                        found = True
                        break
            if found and seek != -1:
                potential_location = seek
            else:
                # if the model let us down, just stick it anywhere
                potential_location = random.randint(region[0], region[1])

        # Record the info for this variant. Listing them one per line to make the indexes later easier to follow
        mutation_data[variant] = [genotype,
                                  region,
                                  is_indel,
                                  is_insertion,
                                  length,
                                  potential_location]
    
    for mutation in mutation_data:
        # Check if we're somewhere we shouldn't be
        allowed = contig_sequence[mutation[5]] in ALLOWED_NUCL

        # Check target and discard bed files
        # Note that if we aren't completely discarding off-target matches, then for the vcf,
        # we don't need to even check if it's in the target region.
        if target_regions and options.discard_offtarget and allowed:
            in_region = False
            for coords in target_regions:
                # Check that this location is valid.
                # Note that we are assuming the coords are half open here.
                if coords[0] <= mutation[5] < coords[1]:
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
                if coords[0] <= mutation[5] < coords[1]:
                    in_region = True
                    break
            if in_region:
                allowed = False

        # Let's assume we're fine, then adjust if not
        final_position = mutation[5]
        # if we're somewhere we shouldn't be, let's find the closest place where we can be.
        if not allowed:
            # See if there is any place after to put it
            sub_sequence = contig_sequence[mutation[5]:]
            # To get the final position relative to the reference, we add subsequence start point
            final_position = min([sub_sequence.index(n) for n in ALLOWED_NUCL])
            if final_position == -1:
                # okay, so we didn't find it to the right of the starting point. To look the other direction
                # We will look up to the current location, but reverse the list, so we find the highest index
                sub_sequence = contig_sequence[:mutation[5]]
                for i in range(len(sub_sequence), 0, 1):
                    if sub_sequence[i] in ALLOWED_NUCL:
                        final_position = i
                # We still didn't find anywhere to put it. I guess we skip it. This seems like an edge case.
                if final_position == -1:
                    print_and_log(f'mutation skipped! {mutation}', 'warning')
                    continue
            # Since we checked a subsequence for the final position, and at this point we know that
            # we have a solid final_position, let's add that to the total.
            final_position += mutation[5]

        # at this point we should have a final position.
        if is_indel:
            if is_insertion:
                ref = contig_sequence[final_position]
                alt = ref + "".join(random.choices(ALLOWED_NUCL, k=length))
            else:
                # plus one because we leave the first base alone in vcf format
                ref = contig_sequence[final_position: length + 1]
                alt = contig_sequence[final_position]
        # if it's not an indel, it's a SNP
        else:
            # See what the ref is
            ref = contig_sequence[final_position]
            # mutate the hell out of it
            alt = random.choice(ALLOWED_NUCL)

        # we're maximally confident that we inserted this variant
        quality_score = max(models.sequencing_error_model.quality_scores)

        # Note that ID and INFO is always just a period, indicating no data, for NEAT
        line = f'{chrom}\t{final_position}\t.\t{ref}\t' \
               f'{alt}\t{quality_score}\t' \
               f'PASS\t.\t' \
               f'GT\t{"/".join(genotype)}\n'

        with open(tmp_vcf_fn, 'a') as tmp:
            tmp.write(line)

        if options.debug:
            print_and_log(f'wrote variant: {line}', 'debug')

    print_and_log(f"Finished mutating {chrom}. Time: {time.time() - start}", 'debug')

    # Let's write out the vcf, if asked to
    if options.produce_vcf:
        print_and_log(f'Writing output vcf', 'info')

        path = f"{tmp_dir_path}/chunk_*.vcf"
        chunksize = 1_000_000
        fid = 1
        lines = []

        with open(str(tmp_vcf_fn), 'r') as f:
            f_out = open(f'{tmp_dir_path}/chunk_{fid}.vcf', 'w')
            for line_num, line in enumerate(f, 1):
                if not line.startswith('@'):
                    lines.append(line)
                if not line_num % chunksize:
                    lines = sorted(lines, key=lambda k: (k.split()[0], int(k.split()[1])))
                    f_out.writelines(lines)
                    print_and_log(f"splitting {chrom} {threadidx} {fid}", 'debug')
                    f_out.close()
                    lines = []
                    fid += 1
                    f_out = open(f'{tmp_dir_path}/chunk_{fid}.vcf' 'w')
            if lines:
                print_and_log(f"splitting {chrom} {threadidx} {fid}", 'debug')
                lines = sorted(lines, key=lambda k: (k.split()[0], int(k.split()[1])))
                f_out.writelines(lines)
                f_out.close()
                lines = []

        chunks = []
        for filename in glob.glob(path):
            chunks += [open(filename, 'r')]

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
