#!/usr/bin/env python

"""
The purpose of this script is to run gen reads in multithreaded mode. The way it will work is to
split the reference fasta up into chromosomes and then run gen_reads in one thread per chromosome mode.
Finally, it will recombine the final output into one.

Note that I'm planninng to hijack NEAT's features for this.
"""

import argparse
import bisect
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
import atexit

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from source.Models import Models
from source.Options import Options
from source.bed_func import parse_bed
from source.ref_func import process_reference
from source.constants_and_defaults import ALLOWED_NUCL
from source.constants_and_defaults import VERSION
from source.error_handling import premature_exit, log_mssg
from source.output_file_writer import OutputFileWriter
from source.run_neat import execute_neat
from source.vcf_func import parse_input_vcf


class StdoutFilter(logging.Filter):
    # Setup filter for logging options
    def filter(self, record):
        return record.levelno == logging.INFO


# useful functions for the rest
def print_start_info():
    """
    Prints a NEAT start message with version number and logs the start time

    :return: logs the time this function was called as start time.
    """
    starttime = time.time()
    print('\n-----------------------------------------------------------')
    log_mssg(f"NEAT multithreaded version, v{VERSION}, is running.", 'info')
    log_mssg(f'Started: {str(starttime)}.', 'info')
    return starttime


def print_end_info(output_class, cancer_output_class, starting_time):
    endtime = time.time()
    list_to_iterate = output_class.files_to_write
    if cancer_output_class:
        list_to_iterate += cancer_output_class.files_to_write
    log_mssg("\nWrote the following files:", 'info')
    i = 1
    for item in list_to_iterate:
        log_mssg(f"\t{i}. {item}", 'info')
        i += 1
    log_mssg("NEAT finished successfully.", "info")
    log_mssg(f"Total runtime: {str(endtime - starting_time)}", 'info')
    print("-------------------------------------------------------------------------")


def print_configuration(args, options):
    """
    Combines the relevant parts of the input args and the options file to print out a
    list of the configuration parameters. Useful for reproducibility.
    """
    log_mssg(f'Run Configuration...', 'INFO')
    potential_filetypes = ['vcf', 'bam', 'fasta', 'fastq']
    log = ''
    for suffix in potential_filetypes:
        key = f'produce_{suffix}'
        if options.args[key]:
            log += f'{suffix} '
    log_mssg(f'Producing the following files: {log.strip()}', 'INFO')
    log_mssg(f'Input file: {options.reference}', 'INFO')
    log_mssg(f'Output files: {args.output}.<{log}>', 'INFO')
    if options.threads == 1:
        log_mssg(f"Single threading - 1 thread.", 'info')
    else:
        log_mssg(f'Multithreading coming soon!!', 'info')
        log_mssg(f"Single threading - 1 thread.", 'info')
        # We'll work on mulitthreading later...
        # print_and_log(f'Multithreading - {options.threads} threads', 'info')
    if options.paired_ended:
        log_mssg(f'Running in paired-ended mode.', 'INFO')
        if options.fragment_model:
            log_mssg(f'Using fragment length model: {options.fragment_model}', 'INFO')
        else:
            log_mssg(f'Using fragment model based on mean={options.fragment_mean}, '
                     f'st dev={options.fragment_st_dev}', 'INFO')
    else:
        log_mssg(f'Running in single-ended mode.', 'INFO')
    log_mssg(f'Using a read length of {options.read_len}', 'INFO')
    log_mssg(f'Average coverage: {options.coverage}', 'INFO')
    log_mssg(f'Using error model: {options.error_model}', 'INFO')
    if options.avg_seq_error:
        log_mssg(f'User defined average sequencing error rate: {options.avg_seq_error}.', 'INFO')
    if options.rescale_qualities:
        log_mssg(f'Quality scores will be rescaled to match avg seq error rate.', 'INFO')
    log_mssg(f'Ploidy value: {options.ploidy}', 'INFO')
    if options.include_vcf:
        log_mssg(f'Vcf of variants to include: {options.include_vcf}', 'INFO')
    if options.target_bed:
        log_mssg(f'BED of regions to target: {options.target_bed}', 'INFO')
        log_mssg(f'Off-target coverage rate: {options.off_target_coverage}', 'INFO')
        log_mssg(f'Discarding off target regions: {options.discard_offtarget}', 'INFO')
    if options.discard_bed:
        log_mssg(f'BED of regions to discard: {options.discard_bed}', 'INFO')
    if options.mutation_model:
        log_mssg(f'Using mutation model in file: {options.mutation_model}', 'INFO')
    if options.mutation_rate:
        log_mssg(f'Rescaling average mutation rate to: {options.mutation_rate}', 'INFO')
    if options.mutation_bed:
        log_mssg(f'BED of mutation rates of different regions: {options.mutation_bed}', 'INFO')
    if options.n_cutoff:
        log_mssg(f'N-cutoff quality score: {options.n_cutoff}', 'INFO')
    if options.gc_model:
        log_mssg(f'Using GC model: {options.gc_model}', 'INFO')
    if options.force_coverage:
        log_mssg(f'Ignoring models and forcing coverage value.', 'INFO')
    log_mssg(f'Debug Mode Activated.', 'debug')
    log_mssg(f'RNG seed value: {options.rng_value}', 'debug')


def find_file_breaks(threads: int, mode: str, reference_index: dict) -> dict:
    """
    Returns a dictionary with the chromosomes as keys
    For the chrom method, the value for each key will just  be "all"
    whereas for subdivison, the value for each key should be a list of indices that partition the
    sequence into roughly equal sizes.
    :param threads: number of threads for this run
    :param mode: partition mode for this run (chrom or subdivision)
    :param reference_index: a dictionary with chromosome keys and sequence values
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

        log_mssg(f'breaks = {partitions}', 'debug')
    else:
        log_mssg("Invalid partition mode. Must be either chrom or subdivision.", 'error')
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


def parse_mutation_rate_dict(mutation_rate_map, avg_rate, non_n_regions):
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
                log_mssg(f'Mutation bed with invalid coordinates', 'error')
                premature_exit(1)
            current_index = regions[i][0]
            # This checks for an overlap. This can only happen when i > 0, by design
            if current_index < prev_index:
                log_mssg(f'Mutation rate region has overlapping regions, merging.', 'warning')
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

    log_file = f'{args.output}.log'
    if pathlib.Path(log_file).exists():
        os.remove(log_file)

    starttime = print_start_info()

    if not args.conf.is_file():
        log_mssg(f'Configuration file ({args.conf}) cannot be found.', 'error')
        premature_exit(1)

    # Reads in the user-entered options from the config file and performs
    # some basic checks and corrections.
    options = Options(args.conf)

    if not options.debug:
        file_handler.setLevel(logging.INFO)

    print_configuration(args, options)

    """
    Model preparation
    
    Read input models or default models, as specified by user.
    """
    log_mssg("Reading Models...", 'info')

    models = Models(options)

    """
    Process Inputs
    """
    log_mssg("Processing inputs...", 'info')

    log_mssg(f'Reading {options.reference}...', 'info')

    reference_index = SeqIO.index(str(options.reference), 'fasta')

    safe_zones, trinucleotide_models = process_reference(reference_index, models)

    log_mssg(f'Reference file indexed.', 'debug')

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
        log_mssg(f"Reading input VCF...", 'info')
        if options.cancer:
            (sample_names, input_variants) = parse_input_vcf(options.include_vcf,
                                                             tumor_normal=True,
                                                             ploidy=options.ploidy)
            # TODO figure out what these were going to be used for
            tumor_ind = sample_names.index('tumor_sample')
            normal_ind = sample_names.index('normal_sample')
        else:
            (sample_names, input_variants) = parse_input_vcf(options.include_vcf,
                                                             ploidy=options.ploidy,
                                                             debug=options.debug)

        log_mssg("Finished reading @include_vcf file.", "info")

    # parse input targeted regions, if present.
    if options.target_bed or options.discard_bed or options.mutation_bed:
        log_mssg(f"Reading input bed files.", 'info')

    # Note if any bed is empty, parse_bed just returns a dict of chromosomes with empty list values

    target_regions_dict = parse_bed(options.target_bed, options.reference_chromosomes,
                                    begins_with_chr, False, options.debug)

    discard_regions_dict = parse_bed(options.discard_bed, options.reference_chromosomes,
                                     begins_with_chr, False, options.debug)

    mutation_rate_dict = parse_bed(options.mutation_bed, options.reference_chromosomes,
                                   begins_with_chr, True, options.debug)

    mutation_rate_dict = parse_mutation_rate_dict(mutation_rate_dict, models.mutation_model['avg_mut_rate'], safe_zones)

    log_mssg(f'Finished reading input beds.', 'debug')

    """
    Initialize Output Files
    """
    log_mssg("Initializing output files...", 'info')

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
        log_mssg(f'Creating output dir: {out_prefix_parent_dir}', 'debug')
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

    log_mssg(f'Output files ready for writing.', 'debug')

    """
    Begin Analysis
    """
    log_mssg("Beginning analysis...", 'info')

    # Find break points in the input file.
    # TODO Debug:
    breaks = find_file_breaks(options.threads, options.partition_mode, reference_index)

    if not breaks:
        # Printing out summary information and end time
        log_mssg("Found no chromosomes in reference.", 'debug')
        print_end_info(output_file_writer, output_file_writer_cancer, starttime)
        sys.exit(0)

    log_mssg("Input file partitioned.", 'debug')

    """
    REDO with mpire. The package can map a function to an iterable, so we will use the breaks as the iterable
    And recast neat as a function (run_neat) that will take in the items from that iterable and process them.
    I think the core of run_neat should be taking in a sequence and returning a mutated sequence. The core of NEAT
    at the moment tries to be everything to everyone.
    """
    log_mssg("Beginning simulation", 'info')
    # these will be the features common to each contig, for multiprocessing
    common_features = {}

    out_prefix = f'{out_prefix_parent_dir}/{out_prefix_name}'

    vcf_files = []

    for contig in breaks:
        inputs_df = pd.DataFrame()
        if not input_variants.empty:
            # We used to filter above, but this has the same effect
            inputs_df = input_variants[input_variants.CHROM.isin([contig])]
        vcf_files.append(execute_neat(reference_index[contig],
                         contig,
                         safe_zones[contig],
                         trinucleotide_models[contig],
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

    now = time.localtime()
    now = time.strftime('%Y_%m_%d_%H%M', now)

    log_dir = os.getcwd()

    log_name = f'{log_dir}/{now}_NEAT.log'

    neat_log = logging.getLogger()
    neat_log.setLevel(logging.DEBUG)
    neat_log.handlers = []

    formatter = logging.Formatter('%(asctime)s %(levelname)s %(threadName)s: %(message)s', datefmt='%Y/%m/%d %H:%M')

    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.addFilter(StdoutFilter())
    stdout_handler.setFormatter(formatter)

    neat_log.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(stream=sys.stderr)
    stderr_handler.setFormatter(formatter)
    stderr_handler.setLevel(logging.WARNING)

    neat_log.addHandler(stderr_handler)

    file_handler = logging.FileHandler(log_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    neat_log.addHandler(file_handler)

    main()
