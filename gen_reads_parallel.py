#!/usr/bin/env python

"""
The purpose of this script is to run gen reads in multithreaded mode. The way it will work is to
split the reference fasta up into chromosomes and then run gen_reads in one thread per chromosome mode.
Finally, it will recombine the final output into one.

Note that I'm planninng to hijack NEAT's features for this.
"""

import os
import sys
import argparse
import datetime
import pathlib
import random
import time
import gzip
import numpy as np
import pickle
import logging
from copy import deepcopy
import multiprocessing
import tempfile
import pandas as pd
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from types import SimpleNamespace
from cpu_load_generator import load_single_core, load_all_cores, from_profile
from mpire import WorkerPool

from source.constants_and_models import VERSION, LOW_COVERAGE_THRESHOLD
from source.Models import parse_input_mutation_model, pickle_load_model
from source.SequencingErrors import SequencingErrors
from source.error_handling import premature_exit, print_and_log

from source.ref_func import find_n_regions
from source.bed_func import parse_bed
from source.vcf_func import parse_vcf
from source.output_file_writer import OutputFileWriter, reverse_complement, sam_flag
from source.probability import DiscreteDistribution, mean_ind_of_weighted_list
from source.SequenceContainer import SequenceContainer
from source.constants_and_models import ALLOWED_NUCL
from source.neat_cigar import CigarString
from source.Models import Models
from source.Tracker import Tracker

# Constants
NEAT_PATH = pathlib.Path(__file__).resolve().parent
WORKING_DIR = pathlib.Path(os.getcwd())


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
    Prints out file names and multithreading info
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
        print_and_log(f'Multithreading - {options.threads} threads', 'info')
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
    if options.debug:
        print_and_log(f'Debug Mode Activated.', 'INFO')
    if options.rng_value:
        print_and_log(f'RNG seed value: {options.rng_value}', 'INFO')


def find_file_breaks(threads: int, mode: str, reference_index: dict, debug: bool) -> dict:
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

    >>> index = {'chr1': "ACCATACAC"}
    >>> find_file_breaks(5, "subdivision", index, False)
    {'chr1': [0, 1, 2, 3, 4]}

    >>> find_file_breaks(5, 'chrom', index, False)
    {'chr1': [0]}
    """
    partitions = {}
    if mode.lower() == "chrom" or threads == 1:
        for contig in reference_index.keys():
            partitions[contig] = [0]
    elif mode.lower() == "subdivision":
        # Instead of this, this should try to subdivide the chrom into chunks
        total_length = 0
        length_dict = {}
        for chrom in reference_index:
            total_length += len(reference_index[chrom])
            length_dict[chrom] = len(reference_index[chrom])

        # Add items one at a time to partition list until the total length is greater than delta.
        index = 0
        for item in length_dict:
            delta = length_dict[item] // threads
            if item not in partitions:
                partitions[item] = []
            for i in range(threads):
                partitions[item].append(index)
                index += delta

        if debug:
            print_and_log(f'breaks = {partitions}', 'debug')
    else:
        print_and_log("Invalid partition mode. Must be either chrom or subdivision.", 'error')
        premature_exit(1)
    return partitions


def find_n_regions(sequence, run_options):
    previous_n_index = 0
    n_count = 0
    n_atlas = []
    for i in range(len(sequence)):
        # Checks for invalid characters, most commonly N
        if sequence[i] not in ALLOWED_NUCL:
            if n_count == 0:
                previous_n_index = i
            n_count += 1
            if i == len(sequence) - 1:
                n_atlas.append((previous_n_index, previous_n_index + n_count))
        else:
            if n_count > 0:
                n_atlas.append((previous_n_index, previous_n_index + n_count))
            n_count = 0


# class representing the options
class Options(SimpleNamespace):
    def __init__(self, config_file):
        self.defs = {}
        self.config_file = config_file

        # Options flags for gen_reads. This metadata dict gives the type of variable (matching the python types)
        # the default value ('.' means no default), and checks. There are four fields: option type (corresponding to
        # a python type), default (or None, if no default), criteria 1 and criteria 2. Any items without values should
        # use None as a placeholder.
        #
        # Option Type: Current possible valuse are string, int, float, and boolean. These correspond to str, int, float,
        # and bool Python types. For files and directories, use string type.
        #
        # Criteria 1: There are two modes of checking: files and numbers. For files, criteria 1 should be set to
        # 'exists' to check file existence or None to skip a check (as with temp_dir, because that is not a user input.
        # For numbers, criteria 1 should be the lowest acceptable value (inclusive) for that variable.
        #
        # Criteria 2: For files, criteria 2 will not be checked, so set to None for consistency. For numbers, this
        # should be the highest acceptable value (inclusive).
        # (type, default, criteria1 (low/'exists'), criteria2 (high/None))
        arbitrarily_large_number = 1e8
        self.defs['reference'] = ('string', None, 'exists', None)
        self.defs['partition_mode'] = ('string', None, None, None)
        self.defs['read_len'] = ('int', 101, 10, arbitrarily_large_number)
        self.defs['threads'] = ('int', 1, 1, arbitrarily_large_number)
        self.defs['coverage'] = ('float', 10.0, 1, arbitrarily_large_number)
        self.defs['error_model'] = ('string', NEAT_PATH / 'models/errorModel_default.pickle.gz', 'exists', None)
        self.defs['avg_seq_error'] = ('float', None, 0, 0.3)
        self.defs['rescale_qualities'] = ('boolean', False, None, None)
        self.defs['ploidy'] = ('int', 2, 1, 100)
        self.defs['include_vcf'] = ('string', None, 'exists', None)
        self.defs['target_bed'] = ('string', None, 'exists', None)
        self.defs['discard_bed'] = ('string', None, 'exists', None)
        self.defs['off_target_coverage'] = ('float', 0.0, 0, 1)
        self.defs['discard_offtarget'] = ('boolean', False, None, None)
        self.defs['mutation_model'] = ('string', None, 'exists', None)
        self.defs['mutation_rate'] = ('float', None, 0, 0.3)
        self.defs['mutation_bed'] = ('string', None, 'exists', None)
        self.defs['n_handling'] = ('string', None, None, None)

        # Params for cancer (not implemented yet)
        self.defs['cancer'] = ('boolean', False, None, None)
        self.defs['cancer_model'] = ('string', None, 'exists', None)
        self.defs['cancer_purity'] = ('float', 0.8, 0.0, 1.0)

        self.defs['n_cutoff'] = ('int', None, 1, 40)
        self.defs['gc_model'] = ('string', NEAT_PATH / 'models/gcBias_default.pickle.gz', 'exists', None)
        self.defs['paired_ended'] = ('boolean', False, None, None)
        self.defs['fragment_model'] = ('string', None, 'exists', None)
        self.defs['fragment_mean'] = ('float', None, 1, arbitrarily_large_number)
        self.defs['fragment_st_dev'] = ('float', None, 1, arbitrarily_large_number)
        self.defs['produce_bam'] = ('boolean', False, None, None)
        self.defs['produce_vcf'] = ('boolean', False, None, None)
        self.defs['produce_fasta'] = ('boolean', False, None, None)
        self.defs['produce_fastq'] = ('boolean', True, None, None)
        self.defs['force_coverage'] = ('boolean', False, None, None)
        self.defs['temp_dir'] = ('string', WORKING_DIR, None, None)
        self.defs['debug'] = ('boolean', False, None, None)
        self.defs['rng_value'] = ('int', None, None, None)

        # Cancer options (not yet implemented)
        self.cancer = False
        self.cancer_model = None
        self.cancer_purity = 0.8

        # Read the config file
        self.args = {}
        self.read()
        # Some options checking to clean up the args dict.
        self.check_options()

        self.__dict__.update(self.args)

    def set_value(self, key, value):
        if key in self.defs.keys():
            self.check_and_log_error(key, value, self.defs[key][2], self.defs[key][3])
        self.args[key] = value
        self.__dict__.update(self.args)

    @staticmethod
    def check_and_log_error(keyname, value_to_check, lowval, highval):
        if lowval != "exists" and highval:
            if not (lowval <= value_to_check <= highval):
                print_and_log(f'@{keyname} must be between {lowval} and {highval}.', 'error')
                premature_exit(1)
        elif lowval == "exists":
            if not pathlib.Path(value_to_check).is_file():
                print_and_log(f'The file given to @{keyname} does not exist', 'error')
                premature_exit(1)
        elif not lowval and not highval:
            # This indicates a boolean or dir and we have nothing to check
            pass
        else:
            print_and_log(f'Problem criteria ({lowval, highval}) in Options definitions for {keyname}.', 'critical')
            premature_exit(1)

    def read(self):
        for line in open(self.config_file):
            line = line.strip()
            if line.startswith('@'):
                line_split = [x.strip().strip('@') for x in line.split('=')]
                # If set to a period, then ignore, by convention.
                if line_split[1] == '.':
                    continue
                key = line_split[0]
                # We can ignore any keys users added but haven't coded for. It must be in the defs dict to be examined.
                if key in list(self.defs.keys()):
                    type_of_var, default, criteria1, criteria2 = self.defs[key]
                    # if it's already set to the default value, ignore.
                    if line_split[1] == default:
                        continue

                    # Now we check that the type is correct and it is in range, depending on the type defined for it
                    # If it passes that it gets put into the args dictionary.
                    if type_of_var == 'string':
                        temp = line_split[1]
                        self.check_and_log_error(key, temp, criteria1, criteria2)
                        self.args[key] = temp
                    elif type_of_var == 'int':
                        try:
                            temp = int(line_split[1])
                        except ValueError:
                            print_and_log(f'The value for {key} must be an integer. No decimals allowed.', 'error')
                            premature_exit(1)
                        self.check_and_log_error(key, temp, criteria1, criteria2)
                        self.args[key] = temp
                    elif type_of_var == 'float':
                        try:
                            temp = float(line_split[1])
                        except ValueError:
                            print_and_log(f'The value for {key} must be a float.', 'error')
                            premature_exit(1)
                        self.check_and_log_error(key, temp, criteria1, criteria2)
                        self.args[key] = temp
                    elif type_of_var == 'boolean':
                        try:
                            if line_split[1].lower() == 'true':
                                self.args[key] = True
                            elif line_split[1].lower() == 'false':
                                self.args[key] = False
                            else:
                                raise ValueError
                        except ValueError:
                            print_and_log(f'\nBoolean key @{key} requires a value of "true" or "false" '
                                          f'(case insensitive).', 'error')
                            premature_exit(1)
                    else:
                        print_and_log(f'BUG: Undefined type in the Options dictionary: {type_of_var}.', 'critical')
                        premature_exit(1)
        # Anything we skipped in the config gets the default value
        # No need to check since these are already CAREFULLY vetted (right!?)
        for key, (_, default, criteria1, criteria2) in self.defs.items():
            if key not in list(self.args.keys()):
                self.args[key] = default

    def check_options(self) -> int:
        """
        Some sanity checks and corrections to the options.
        """

        if self.args['produce_fasta']:
            print_and_log("\nFASTA mode active.", 'info')
            print_and_log("NOTE: At the moment, NEAT can produce a FASTA or FASTQ files, not both.", 'info')
            # Turn off fastq and paired-ended mode for now
            self.args['produce_fastq'] = False
            self.args['paired_ended'] = False
            self.args['fragment_model'] = None
            self.args['fragment_mean'] = None
            self.args['fragment_st_dev'] = None
        if not self.args['produce_bam'] and not self.args['produce_vcf'] \
                and not self.args['produce_fasta'] and not self.args['produce_fastq']:
            print_and_log('No files would be produced, as all file types are set to false', 'error')
            premature_exit(1)
        if not self.args['produce_fastq']:
            print_and_log("Bypassing FASTQ generation.", 'info')
        if self.args['produce_vcf'] and (not self.args['produce_fastq'] and not self.args['produce_bam']
                                         and not self.args['produce_fasta']):
            print_and_log('Only producing VCF output.', 'info')
        if self.args['produce_bam'] and (not self.args['produce_fastq'] and not self.args['produce_vcf']
                                         and not self.args['produce_fasta']):
            print_and_log('Only producing BAM output.', 'info')

        # This next section just checks all the paired ended stuff
        flagged = False
        if self.args['paired_ended']:
            if self.args['fragment_model']:
                self.args['fragment_mean'] = None
                self.args['fragment_st_dev'] = None
            elif self.args['fragment_mean']:
                if not self.args['fragment_st_dev']:
                    flagged = True
            else:
                flagged = True
        if flagged:
            print_and_log("For paired ended mode, you need either a "
                          "@fragment_model or both @fragment_mean and @fragment_st_dev", 'error')
            premature_exit(1)

        self.args['n_handling'] = "ignore"
        if self.args['paired_ended']:
            self.args['n_handling'] = "random"

# command line interface
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

    # Set the random seed. If rng_value is None, then the default for random.seed is to use system time.
    random.seed(options.rng_value)

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
    # TODO still not sure I need this. It requires combing through the entire sequence, which I don't like
    # Maybe add this into the main loop.
    # n_regions = find_n_regions(reference_index, options)

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

    input_variants = None
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

    # parse input targeted regions, if present
    if options.target_bed or options.discard_bed or options.mutation_bed:
        print_and_log(f"Reading input bed files.", 'info')
    # Note if options.target_bed is empty, this just returns an empty data frame
    target_regions_df = parse_bed(options.target_bed, options.reference_chromosomes,
                                  begins_with_chr, False, options.debug)

    # parse discard bed similarly
    discard_regions_df = parse_bed(options.discard_bed, options.reference_chromosomes,
                                   begins_with_chr, False, options.debug)

    # parse input mutation rate rescaling regions, if present
    mutation_rate_df = parse_bed(options.mutation_bed, options.reference_chromosomes,
                                 begins_with_chr, True, options.debug)

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

    # Initialize simulation
    thread_index = 0
    processes = []

    """
    REDO with mpire. The package can map a function to an iterable, so we will use the breaks as the iterable
    And recast neat as a function (run_neat) that will take in the items from that iterable and process them.
    I think the core of run_neat should be taking in a sequence and returning a mutated sequence. Currently, NEAT
    """
    for chrom in breaks:
        thread_index += 1


    # Step 4 (CAVA 496 - 497) - Merging tmp files and writing out final files


    # Step 5 (CAVA 500 - 501) - Printing out summary and end time


if __name__ == '__main__':
    main()
