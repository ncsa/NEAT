#!/usr/bin/env python

"""
The purpose of this script is to run gen reads in multithreaded mode. The way it will work is to
split the reference fasta up into chromosomes and then run gen_reads in one thread per chromosome mode.
Finally, it will recombine the final output into one.

Note that I'm planninng to hijack NEAT's features for this.
"""

import argparse
import bisect
import datetime
import gzip
import logging
import pathlib
import random
import sys
import time
import pandas as pd
import glob
from heapq import merge
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from source.Models import Models
from source.Options import Options
from source.bed_func import parse_bed
from source.constants_and_defaults import ALLOWED_NUCL
from source.constants_and_defaults import VERSION
from source.error_handling import premature_exit, print_and_log
from source.output_file_writer import OutputFileWriter
from source.run_neat import execute_neat
from source.vcf_func import parse_vcf


# useful functions for the rest
def print_start_info() -> datetime.datetime:
    """
    Prints a NEAT start message with version number and logs the start time

    :return: logs the time this function was called as start time.
    """
    starttime = datetime.datetime.now()
    print('\n-----------------------------------------------------------')
    print_and_log(f"NEAT multithreaded version, v{VERSION}, is running.", 'info')
    print_and_log(f'Started: {str(starttime)}.', 'info')
    return starttime


def print_end_info(output_class, cancer_output_class, starting_time):
    endtime = datetime.datetime.now()
    list_to_iterate = output_class.files_to_write
    if cancer_output_class:
        list_to_iterate += cancer_output_class.files_to_write
    print_and_log("\nWrote the following files:", 'info')
    i = 1
    for item in list_to_iterate:
        print_and_log(f"\t{i}. {item}", 'info')
        i += 1
    print_and_log("NEAT finished successfully.", "info")
    print_and_log(f"Total runtime: {str(endtime-starting_time)}", 'info')
    print("-------------------------------------------------------------------------")


def print_configuration(args, options):
    """
    Combines the relevant parts of the input args and the options file to print out a
    list of the configuration parameters. Useful for reproducibility.
    """
    print_and_log(f'Run Configuration...', 'INFO')
    potential_filetypes = ['vcf', 'bam', 'fasta', 'fastq']
    log = ''
    for suffix in potential_filetypes:
        key = f'produce_{suffix}'
        if options.args[key]:
            log += f'{suffix} '
    print_and_log(f'Producing the following files: {log.strip()}', 'INFO')
    print_and_log(f'Input file: {options.reference}', 'INFO')
    print_and_log(f'Output files: {args.output}.<{log}>', 'INFO')
    if options.threads == 1:
        print_and_log(f"Single threading - 1 thread.", 'info')
    else:
        print_and_log(f'Multithreading coming soon!!', 'info')
        premature_exit(0)
        # We'll work on mulitthreading later...
        # print_and_log(f'Multithreading - {options.threads} threads', 'info')
    if options.paired_ended:
        print_and_log(f'Running in paired-ended mode.', 'INFO')
        if options.fragment_model:
            print_and_log(f'Using fragment length model: {options.fragment_model}', 'INFO')
        else:
            print_and_log(f'Using fragment model based on mean={options.fragment_mean}, '
                          f'st dev={options.fragment_st_dev}', 'INFO')
    else:
        print_and_log(f'Running in single-ended mode.', 'INFO')
    print_and_log(f'Using a read length of {options.read_len}', 'INFO')
    print_and_log(f'Average coverage: {options.coverage}', 'INFO')
    print_and_log(f'Using error model: {options.error_model}', 'INFO')
    if options.avg_seq_error:
        print_and_log(f'User defined average sequencing error rate: {options.avg_seq_error}.', 'INFO')
    if options.rescale_qualities:
        print_and_log(f'Quality scores will be rescaled to match avg seq error rate.', 'INFO')
    print_and_log(f'Ploidy value: {options.ploidy}', 'INFO')
    if options.include_vcf:
        print_and_log(f'Vcf of variants to include: {options.include_vcf}', 'INFO')
    if options.target_bed:
        print_and_log(f'BED of regions to target: {options.target_bed}', 'INFO')
        print_and_log(f'Off-target coverage rate: {options.off_target_coverage}', 'INFO')
        print_and_log(f'Discarding off target regions: {options.discard_offtarget}', 'INFO')
    if options.discard_bed:
        print_and_log(f'BED of regions to discard: {options.discard_bed}', 'INFO')
    if options.mutation_model:
        print_and_log(f'Using mutation model in file: {options.mutation_model}', 'INFO')
    if options.mutation_rate:
        print_and_log(f'Rescaling average mutation rate to: {options.mutation_rate}', 'INFO')
    if options.mutation_bed:
        print_and_log(f'BED of mutation rates of different regions: {options.mutation_bed}', 'INFO')
    if options.n_cutoff:
        print_and_log(f'N-cutoff quality score: {options.n_cutoff}', 'INFO')
    if options.gc_model:
        print_and_log(f'Using GC model: {options.gc_model}', 'INFO')
    if options.force_coverage:
        print_and_log(f'Ignoring models and forcing coverage value.', 'INFO')
    print_and_log(f'Creating temp files in directory: {options.temp_dir}', 'info')
    if options.debug:
        print_and_log(f'Debug Mode Activated.', 'INFO')
    if options.rng_value:
        print_and_log(f'RNG seed value: {options.rng_value}', 'INFO')


def find_file_breaks(threads: int, mode: str, reference_index: dict, debug: bool = False) -> dict:
    """
    Returns a dictionary with the chromosomes as keys
    For the chrom method, the value for each key will just  be "all"
    whereas for subdivison, the value for each key should be a list of indices that partition the
    sequence into roughly equal sizes.
    :param threads: number of threads for this run
    :param mode: partition mode for this run (chrom or subdivision)
    :param reference_index: a dictionary with chromosome keys and sequence values
    :param debug: Turns debug mode on, if True
    :return: a dictionary containing the chromosomes as keys and either "all" for valuse, or a list of indices

    >>> index = {'chr1': "ACCATACACGGGCAACACACGTACACATTATACC"}
    >>> find_file_breaks(5, "subdivision", index, False)
    {'chr1': [(0, 6), (6, 12), (12, 18), (18, 24), (24, 34)]}
    >>> find_file_breaks(4, 'subdivision', index, False)
    {'chr1': [(0, 8), (8, 16), (16, 24), (24, 34)]}
    >>> find_file_breaks(5, 'chrom', index, False)
    {'chr1': [(0, 34)]}
    >>> find_file_breaks(2, 'subdivision', index, False)
    {'chr1': [(0, 17), (17, 34)]}
    """
    partitions = {}
    if mode.lower() == "chrom" or threads == 1:
        for contig in reference_index.keys():
            partitions[contig] = [(0, len(reference_index[contig]))]
    elif mode.lower() == "subdivision":
        # Add items one at a time to partition list until the total length is greater than delta.
        for contig in reference_index:
            contig_length = len(reference_index[contig])
            delta = contig_length // threads

            if contig not in partitions:
                partitions[contig] = []

            breakpoints = list(range(0, contig_length, delta))

            # Since we know the first index will be zero, we can skip the first item in the breakpoints list
            # And since we want the last partition to grab the rest, we'll stop short and add it manually
            for index in breakpoints:
                if index + delta <= contig_length:
                    partitions[contig].append((index, index + delta))
                else:
                    # Have to extend the last one so we don't get a tiny read we can't process
                    partitions[contig][-1] = (partitions[contig][-1][0], contig_length)

        if debug:
            print_and_log(f'breaks = {partitions}', 'debug')
    else:
        print_and_log("Invalid partition mode. Must be either chrom or subdivision.", 'error')
        premature_exit(1)
    return partitions


def map_reference(index):
    """
    Makes a dictionary with n regions, not n regions and the total n count
    Note that for this program, anything not ACTG is an N
    :param index: Input reference_index biopython object
    :return: A dictionry indexed by contig in the input index, with a subdicitonary containing
             'n_regions' (list of tuples) - start: end for regions containing unallowed characters in contig
             'non_n' (list of tuples) - start: end for regions containing allowed characters in contig
             'n_count' (int) - number of unallowed characters in contig

    >>> my_seq = {'FAKE1': SeqRecord(Seq("ANNUNACTGNNNGTCA"))}
    >>> print(map_reference(my_seq))
    {'FAKE1': {'n_regions': [(1, 5), (9, 12)], 'non_n': [(0, 1), (5, 9)], 'n_count': 7}}
    >>> my_seq = {'FAKE1': SeqRecord(Seq("NNUNACTGNNNGTCA"))}
    >>> print(map_reference(my_seq))
    {'FAKE1': {'n_regions': [(0, 4), (8, 11)], 'non_n': [(4, 8)], 'n_count': 7}}
    """
    reference_atlas = {x: {'n_regions': [], 'non_n': [], 'n_count': 0} for x in list(index.keys())}
    prev_non_n = 0
    prev_index = 0
    n_count = 0
    for contig in index:
        contig_sequence = index[contig].seq
        # Check if the sequence startswith an unallowed character
        startswith = any([contig_sequence.startswith(n) for n in ALLOWED_NUCL])
        for i in range(len(contig_sequence)):
            # Anything not in the allowed list gets called an N
            if contig_sequence[i] not in ALLOWED_NUCL:
                reference_atlas[contig]['n_count'] += 1
                # In the case where the input sequence has Ns in the middle (AAAANNNAAAA),
                # startswith should catch and put the first interval in (0, 4)
                if prev_non_n > 0 or startswith:
                    reference_atlas[contig]['non_n'].append((prev_non_n, i))
                    prev_non_n = 0
                    startswith = False
                if n_count == 0:
                    prev_index = i
                n_count += 1
                if i == len(contig_sequence) - 1:
                    reference_atlas[contig]['n_regions'].append((prev_index, prev_index + n_count))
            else:
                if n_count > 0:
                    reference_atlas[contig]['n_regions'].append((prev_index, prev_index + n_count))
                    prev_non_n = prev_index + n_count
                n_count = 0

    return reference_atlas


def parse_mutation_rate_dict(mutation_rate_map, avg_rate):
    """
    This parses the mutation rate dict, in order to fill in the dict so it can by more easily cycled through
    later.

    TODO write tests:
    - intervals that overlap
    - intervals that don't overlap
    - intervals with gaps
    - intervals with no gaps
    - mixes of the above
    - empty dictionary ({'chrX': []})
    """
    # creates the default dict
    ret_dict = {x: [] for x in mutation_rate_map}
    for chrom in mutation_rate_map:
        regions = mutation_rate_map[chrom]
        regions.sort()
        total_regions = []
        # check for overlaps
        prev_index = 0
        for i in range(len(regions)):
            # This compares the set of coordinates
            if regions[i][0] >= regions[i][1]:
                print_and_log(f'Mutation bed with invalid coordinates', 'error')
                premature_exit(1)
            current_index = regions[i][0]
            # This checks for an overlap. This can only happen when i > 0, by design
            if current_index < prev_index:
                print_and_log(f'Mutation rate region has overlapping regions, merging.', 'warning')
                # We just average the two regions together
                interval = (regions[i - 1][0], regions[i][1], (regions[i - 1][2] + regions[i][2]) / 2)
                # Get rid of previous, to be replaced by ^^
                total_regions.pop()
            else:
                interval = (current_index, regions[i][1], regions[i][2])
                # This creates an interval in between the bed intervals with default mutation rate
                if current_index != prev_index:
                    total_regions.append((prev_index, current_index, avg_rate))
            total_regions.append(interval)
            prev_index = regions[i][1]
        ret_dict[chrom] = total_regions
    return ret_dict


def main(raw_args=None):
    """
    Parallel main function. Takes args and parses the ones needed to start multiprocessing. The rest will get passed
    along to the main processing part of gen_reads.
    """

    parser = argparse.ArgumentParser(description=f'This script runs gen_reads v{VERSION}',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-c', '--config', default='neat_config.txt', dest='conf', metavar='config',
                        help="Any extra arguments valid for NEAT")
    parser.add_argument('-o', '--output', required=True, dest='output', metavar='output',
                        help="Prefix for the output. Can include a path before it.")

    args = parser.parse_args(raw_args)

    # This will by default look for a file called "neat_config.txt" in the current working dir
    args.conf = pathlib.Path(args.conf)

    log_name = f'{args.output}.log'

    logging.basicConfig(filename=log_name, filemode='w',
                        format='%(asctime)s %(levelname)s %(threadName)s: %(message)s',
                        level=logging.DEBUG)

    starttime = print_start_info()

    if not args.conf.is_file():
        print_and_log(f'Configuration file ({args.conf}) cannot be found.', 'error')
        premature_exit(1)

    # Reads in the user-entered options from the config file and performs
    # some basic checks and corrections.
    options = Options(args.conf)
    print_configuration(args, options)

    """
    Model preparation
    
    Read input models or default models, as specified by user.
    """
    print_and_log("Reading Models...", 'info')

    models = Models(options)

    """
    Process Inputs
    """
    print_and_log("Processing inputs...", 'info')

    print_and_log(f'Reading {options.reference}...', 'info')

    reference_index = SeqIO.index(str(options.reference), 'fasta')

    if options.debug:
        print_and_log(f'Reference file indexed.', 'debug')

    # there is not a reference_chromosome in the options defs because
    # there is nothing that really needs checking, so this command will just
    # set the value. This might be an input in the future, for example if you want
    # to only simulate certain chroms.
    options.set_value('reference_chromosomes', list(reference_index.keys()))

    # TODO maybe take this out. I feel like we're overthinking chr, but I could be wrong
    # this parameter will ensure we account for 'chr' in the output.
    # It is True if they all start with "chr" and False otherwise.
    begins_with_chr = all(k.startswith('chr') for k in options.reference_chromosomes)

    sample_names = []
    input_variants = pd.DataFrame()
    if options.include_vcf:
        print_and_log(f"Reading input VCF...", 'info')
        if options.cancer:
            (sample_names, input_variants) = parse_vcf(options.include_vcf,
                                                       tumor_normal=True,
                                                       ploidy=options.ploidy,
                                                       debug=options.debug)
            # TODO figure out what these were going to be used for
            tumor_ind = sample_names.index('tumor_sample_split')
            normal_ind = sample_names.index('normal_sample_split')
        else:
            (sample_names, input_variants) = parse_vcf(options.include_vcf,
                                                       ploidy=options.ploidy,
                                                       debug=options.debug)

        if options.debug:
            print_and_log("Finished reading @include_vcf file. Now filtering.", "debug")

        # Remove any chromosomes that aren't in the reference.
        input_variants_chroms = input_variants['CHROM'].unique()
        for item in input_variants_chroms:
            if item not in options.reference_chromosomes and not printed_warning:
                print_and_log(f'Warning: ignoring all input vcf records for {item} '
                              f'because it is not found in the reference.', 'warning')
                print_and_log(f'\tIf this is unexpected, check that that {item} '
                              f'matches reference name exactly.', 'warning')
                printed_warning = True
                input_variants = input_variants[input_variants['CHROM'] != item]

        # Check the variants and classify as needed
        for chrom in options.reference_chromosomes:
            n_skipped = [0, 0, 0]
            if chrom in input_variants_chroms:
                for index, row in input_variants[input_variants['CHROM'] == chrom].iterrows():
                    span = (row['POS'], row['POS'] + len(row['REF']))
                    # -1 because going from VCF coords to array coords
                    r_seq = str(reference_index[chrom].seq[span[0] - 1:span[1] - 1])
                    # Checks if there are any invalid nucleotides in the vcf items
                    any_bad_nucl = any((nn not in ALLOWED_NUCL) for nn in
                                       [item for sublist in row['alt_split'] for item in sublist])
                    # Ensure reference sequence matches the nucleotide in the vcf
                    if r_seq != row['REF']:
                        n_skipped[0] += 1
                        input_variants.drop(index, inplace=True)
                        continue
                    # Ensure that we aren't trying to insert into an N region
                    elif 'N' in r_seq:
                        n_skipped[1] += 1
                        input_variants.drop(index, inplace=True)
                        continue
                    # Ensure that we don't insert any disallowed characters
                    elif any_bad_nucl:
                        n_skipped[2] += 1
                        input_variants.drop(index, inplace=True)
                        continue

                if options.debug:
                    print_and_log("Finished filtering @include_vcf file.", 'debug')

                print_and_log(f'Found {len(input_variants)} valid variants after filtering for {chrom} in @include_vcf.', 'info')
                if any(n_skipped):
                    print_and_log(f'variants skipped: {sum(n_skipped)}', 'info')
                    print_and_log(f' - [{str(n_skipped[0])}] ref allele did not match reference', 'info')
                    print_and_log(f' - [{str(n_skipped[1])}] attempted to insert into N-region', 'info')
                    print_and_log(f' - [{str(n_skipped[2])}] alt allele contained non-ACGT characters', 'info')
                else:
                    print_and_log(f'variants skipped: 0', 'info')

    # parse input targeted regions, if present.
    if options.target_bed or options.discard_bed or options.mutation_bed:
        print_and_log(f"Reading input bed files.", 'info')

    # Note if any bed is empty, parse_bed just returns a dict of chromosomes with empty list values

    target_regions_dict = parse_bed(options.target_bed, options.reference_chromosomes,
                                    begins_with_chr, False, options.debug)

    discard_regions_dict = parse_bed(options.discard_bed, options.reference_chromosomes,
                                     begins_with_chr, False, options.debug)

    mutation_rate_dict = parse_bed(options.mutation_bed, options.reference_chromosomes,
                                   begins_with_chr, True, options.debug)

    mutation_rate_dict = parse_mutation_rate_dict(mutation_rate_dict, models.mutation_model['avg_mut_rate'])

    if options.debug:
        print_and_log(f'Finished reading input beds.', 'debug')

    """
    Initialize Output Files
    """
    print_and_log("Initializing output files...", 'info')

    # Prepare headers
    bam_header = None
    vcf_header = None
    if options.produce_bam:
        bam_header = reference_index
    if options.produce_vcf:
        vcf_header = [options.reference]

    out_prefix_name = pathlib.Path(args.output).resolve().name
    out_prefix_parent_dir = pathlib.Path(pathlib.Path(args.output).resolve().parent)
    if not out_prefix_parent_dir.is_dir():
        if options.debug:
            print_and_log(f'Creating output dir: {out_prefix_parent_dir}', 'info')
        out_prefix_parent_dir.mkdir(parents=True, exist_ok=True)

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf.
    # We'll also keep track here of what files we are producing.
    if options.cancer:
        output_normal = out_prefix_parent_dir / f'{out_prefix_name}_normal'
        output_tumor = out_prefix_parent_dir / f'{out_prefix_name}_tumor'
        output_file_writer = OutputFileWriter(output_normal,
                                              bam_header=bam_header,
                                              vcf_header=vcf_header,
                                              options_file=options
                                              )
        output_file_writer_cancer = OutputFileWriter(output_tumor,
                                                     bam_header=bam_header,
                                                     vcf_header=vcf_header,
                                                     options_file=options
                                                     )
    else:
        outfile = out_prefix_parent_dir / out_prefix_name
        output_file_writer = OutputFileWriter(outfile,
                                              bam_header=bam_header,
                                              vcf_header=vcf_header,
                                              options_file=options
                                              )
        output_file_writer_cancer = None

    if options.debug:
        print_and_log(f'Output files ready for writing.', 'debug')

    """
    Begin Analysis
    """
    print_and_log("Beginning analysis...", 'info')

    # Find break points in the input file.
    # TODO Debug:
    breaks = find_file_breaks(options.threads, options.partition_mode, reference_index, options.debug)

    if not breaks:
        # Printing out summary information and end time
        if options.debug:
            print_and_log("Found no chromosomes in reference.", 'debug')
        print_end_info(output_file_writer, output_file_writer_cancer, starttime)
        sys.exit(0)

    if options.debug:
        print_and_log("Input file partitioned.", 'debug')

    """
    REDO with mpire. The package can map a function to an iterable, so we will use the breaks as the iterable
    And recast neat as a function (run_neat) that will take in the items from that iterable and process them.
    I think the core of run_neat should be taking in a sequence and returning a mutated sequence. The core of NEAT
    at the moment tries to be everything to everyone.
    """
    print_and_log("Beginning simulation", 'info')
    # these will be the features common to each contig, for multiprocessing
    common_features = {}

    out_prefix = f'{out_prefix_parent_dir}/{out_prefix_name}'

    if options.produce_vcf:
        print_and_log(f"Generating list of mutations (golden vcf)", "info")

    vcf_files = []

    for contig in breaks:
        inputs_df = pd.DataFrame()
        if not input_variants.empty:
            inputs_df = input_variants[input_variants.CHROM.isin([chrom])]
        vcf_files.append(execute_neat(reference_index[contig],
                         contig,
                         out_prefix_name,
                         target_regions_dict[contig],
                         discard_regions_dict[contig],
                         mutation_rate_dict[contig],
                         inputs_df,
                         models,
                         options,
                         out_prefix))

    if options.produce_vcf:
        chunks = []
        # TODO double check that breaks is in the same order as the input fasta
        for vcf in vcf_files:
            chunks += [gzip.open(vcf, 'r')]

        with gzip.open(output_file_writer.vcf_fn, 'ab') as vcf_out:
            vcf_out.writelines(merge(*chunks))

        for item in chunks:
            item.close()
            os.remove(item.name)


    # Step 4 (CAVA 496 - 497) - Merging tmp files and writing out final files

    # Step 5 (CAVA 500 - 501) - Printing out summary and end time

    # End info
    print_end_info(output_file_writer, output_file_writer_cancer, starttime)


if __name__ == '__main__':
    main()
