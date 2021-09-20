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
import gzip
import numpy as np
import pickle
import logging
from copy import deepcopy
import multiprocessing
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from types import SimpleNamespace

from source.constants_and_models import VERSION, LOW_COVERAGE_THRESHOLD
from source.input_file_reader import parse_input_mutation_model
from source.ReadContainer import ReadContainer
from source.error_handling import will_exit, print_and_log

from source.ref_func import find_n_regions
from source.bed_func import parse_bed
from source.vcf_func import parse_vcf
from source.output_file_writer import OutputFileWriter, reverse_complement, sam_flag
from source.probability import DiscreteDistribution, mean_ind_of_weighted_list
from source.SequenceContainer import SequenceContainer
from source.constants_and_models import ALLOWED_NUCL
from source.neat_cigar import CigarString

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
    print(f'INFO - Started: {str(starttime)}.')
    return starttime


def run(opts):
    """

    """
    pass


def find_chromosomes(reference: pathlib.Path):
    """

    """
    pass


def read_chroms(reference: pathlib.Path):
    """
    Input is a fasta reference file. The output is a list of chromosomes.
    """
    pass


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
            log += f' {suffix} '
    print_and_log(f'Producing the following files: {log.strip()}', 'INFO')
    print_and_log(f'Input file: {options.reference}', 'INFO')
    print_and_log(f'Output files: {args.output}.<{log}>', 'INFO')
    if options.threads == 1:
        print_and_log(f"Single threading - 1 thread.", 'info')
    else:
        print_and_log(f'Multithreading - {options.threads} threads', 'info')
    logging.info(f'Threads: {options.threads}')
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


def pickle_load_model(file, mssg) -> list:
    try:
        return pickle.load(file)
    except IOError as e:
        print_and_log(mssg, 'error')
        logging.error(str(e))
        will_exit(1)
    except EOFError as e:
        print_and_log(mssg, 'error')
        logging.error(str(e))
        will_exit(1)
    except ValueError as e:
        print_and_log(mssg, 'error')
        logging.error(str(e))
        will_exit(1)


# class representing the options
class Options(SimpleNamespace):
    def __init__(self, config_file):
        self.defs = {}
        self.config_file = config_file

        # Options flags for gen_reads. This metadata dict gives the type of variable (matching the python types)
        # the default value ('.' means no default), and checks. Files only have one criteria: 'exists' or
        # nothing. If there is nothing after the default value, the code will not check the value at all.
        # If the criteria is 'exists' it will treat it as a file and check that it exists with pathlib.
        # For numbers you need two criteria: a low value and a high value for the range that the variable
        # should fall between. If there is no default and/or no criteria, use None as a placeholder.
        # (type, default, criteria1 (low), criteria2 (high))
        arbitrarily_large_number = 1000000
        self.defs['reference'] = ('string', None, 'exists', None)
        self.defs['read_len'] = ('int', 101, 10, arbitrarily_large_number)
        self.defs['threads'] = ('int', 1, 1, arbitrarily_large_number)
        self.defs['coverage'] = ('float', 10.0, 1, arbitrarily_large_number)
        self.defs['error_model'] = ('string', NEAT_PATH / 'models/errorModel_default.p', 'exists', None)
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
        self.defs['n_cutoff'] = ('int', None, 1, 40)
        self.defs['gc_model'] = ('string', NEAT_PATH / 'models/gcBias_default.p', 'exists', None)
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

        # Read the config file
        self.args = {}
        self.read()
        # Some options checking to clean up the args dict
        self.check_options()

        self.__dict__.update(self.args)

    @staticmethod
    def check_and_log_error(keyname, value_to_check, lowval, highval):
        if lowval != "exists" and highval:
            if not (lowval <= value_to_check <= highval):
                print_and_log(f'@{keyname} must be between {lowval} and {highval}.', 'error')
                will_exit(1)
        elif lowval == "exists":
            if not pathlib.Path(value_to_check).is_file():
                print_and_log(f'the file given to @{keyname} does not exist', 'error')
                will_exit(1)
        elif not lowval and not highval:
            # This indicates a boolean and we have nothing to check
            pass
        else:
            print_and_log(f'Undeclared criteria {lowval, highval} in Options definitions.', 'error')
            will_exit(1)

    def read(self):
        for line in open(self.config_file):
            line = line.strip()
            if line.startswith('@'):
                line_split = [x.strip().strip('@') for x in line.split('=')]
                # If set to a period, then ignore, by convention.
                if line_split[1] == '.':
                    continue
                key = line_split[0]
                # We can ignore any keys users. added but haven't coded for.
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
                            will_exit(1)
                        self.check_and_log_error(key, temp, criteria1, criteria2)
                        self.args[key] = temp
                    elif type_of_var == 'float':
                        try:
                            temp = float(line_split[1])
                        except ValueError:
                            print_and_log(f'The value for {key} must be a float.', 'error')
                            will_exit(1)
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
                            will_exit(1)
                    else:
                        print_and_log(f'BUG: Undefined type in the Options dictionary: {type_of_var}.', 'critical')
                        will_exit(1)
        # Anything we skipped in the config gets the default value
        # No need to check since these are already CAREFULLY vetted
        for key, (_, default, criteria1, criteria2) in self.defs.items():
            if key not in list(self.args.keys()):
                # Let's double check that the file structure is as expected
                if default:
                    self.check_and_log_error(key, default, criteria1, criteria2)
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
            will_exit(1)
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
            print_and_log("\nPaired-ended mode", 'info')
            if self.args['fragment_model']:
                print_and_log(f"Using fragment length model {self.args['fragment_model']} to produce paired ended reads")
                self.args['fragment_mean'] = None
                self.args['fragment_st_dev'] = None
            elif self.args['fragment_mean']:
                if self.args['fragment_st_dev']:
                    print_and_log(f"Using fragment length mean = {self.args['fragment_mean']}, "
                                  f"standard deviation = {self.args['fragment_st_dev']} "
                                  f"to produce paired ended reads.", 'info')
                else:
                    flagged = True
            else:
                flagged = True
        if flagged:
            print_and_log("For paired ended mode, you need either a "
                          "@fragment_model or both @fragment_mean and @fragment_st_dev", 'error')
            will_exit(1)


# class representing a chromosome mutation
class SingleJob(multiprocessing.Process):
    def __init__(self, threadidx, options, copts, endline, genelist, transcriptlist, snplist,
                 impactdir, num_of_records):
        multiprocessing.Process.__init__(self)

        # TODO check into multiprocessing.Pool()


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
                        format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.DEBUG)

    starttime = print_start_info()

    if not args.conf.is_file():
        print_and_log(f'Configuration file ({args.conf}) cannot be found.', 'error')
        will_exit(1)

    # Reads in the user-entered options from the config file and performs
    # some basic checks and corrections.
    options = Options(args.conf)
    print_configuration(args, options)

    # Set the random seed. If rng_value is None, then the default for random.seed is to use system time.
    random.seed(options.rng_value)

    # cancer params (disabled currently)
    options.cancer = False
    options.cancer_model = None
    options.cancer_purity = 0.8

    """
    Model preparation
    """
    if options.debug:
        logging.info("Model preparation")
    # load mutation models:
    mutation_model = parse_input_mutation_model(options.mutation_model, 1)
    if options.cancer:
        cancer_model = parse_input_mutation_model(options.cancer_model, 2)

    # Implements sequencing error model
    mssg = "Problem loading Sequencing error model data @error_model"
    error_data = pickle_load_model(open(options.error_model, 'rb'), mssg)
    sequencing_error_class = ReadContainer(options.read_len, options.error_model,
                                           options.avg_seq_error)
    if options.debug:
        logging.info(f'Sequencing read container created')

    # GC bias model
    mssg = "f'ERROR: problem reading @gc_model. Please check file path and try again. " \
           "This file should be the output of compute_gc.py'"
    gc_scale_count, gc_scale_value = pickle_load_model(open(options.gc_model, 'rb'), mssg)

    gc_window_size = gc_scale_count[-1]

    # Handle paired-ended data, if applicable
    if options.paired_ended:
        fraglen_distribution = None
        if options.fraglen_model:
            mssg = f'ERROR: Problem loading the empirical fragment length model.'
            print("Using empirical fragment length distribution")
            logging.info("Using empirical fragment length distribution")
            potential_values, potential_prob = pickle_load_model(open(options.fraglen_model, 'rb'), mssg)

            fraglen_values = []
            fraglen_probability = []
            for i in range(len(potential_values)):
                if potential_values[1] > options.read_len:
                    fraglen_values.append(potential_values[i])
                    fraglen_probability.append(potential_prob[i])

            fraglen_distribution = DiscreteDistribution(fraglen_probability, fraglen_values)
            fragment_mean = fraglen_values[mean_ind_of_weighted_list(fraglen_probability)]

        # Using artificial fragment length distribution, if the parameters were specified
        # fragment length distribution: normal distribution that goes out to +- 6 standard deviations
        else:
            print(f'Using artificial fragment length distribution.')
            logging.info(f'Using artificial fragment length distribution.')
            if options.fragment_st_dev == 0:
                fraglen_distribution = DiscreteDistribution([1], [options.fragment_mean],
                                                            degenerate_val=options.fragment_mean)
            else:
                potential_values = range(options.fragment_mean - 6 * options.fragment_st_dev,
                                         options.fragment_mean + 6 * options.fragment_st_dev + 1)
                fraglen_values = []
                for i in range(len(potential_values)):
                    if potential_values[i] > options.read_len:
                        fraglen_values.append(potential_values[i])
                fraglen_probability = [np.exp(-(((n - options.fragment_mean) ** 2) /
                                                (2 * (options.fragment_st_dev ** 2)))) for n in
                                       fraglen_values]
                fraglen_distribution = DiscreteDistribution(fraglen_probability, fraglen_values)
        if options.debug:
            logging.info(f'Processed paired-end models')

    """
    Process Inputs
    """
    if options.debug:
        logging.info("Process Inputs")

    print_and_log(f'Reading {options.reference}...', 'info')

    ref_index = SeqIO.index(str(options.reference), 'fasta')
    reference_chromosomes = list(ref_index.keys())
    # this parameter will ensure we account for 'chr' in the output.
    # It is True if they all start with "chr" and False otherwise.
    begins_with_chr = all(k.startswith('chr') for k in reference_chromosomes)

    if options.debug:
        logging.info(f'Reference file indexed.')

    if options.paired_ended:
        n_handling = ('random', options.fragment_mean)
    else:
        n_handling = ('ignore', options.read_len)

    input_variants = None
    printed_warning = False
    if options.include_vcf:
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
            logging.info("Finished reading @include_vcf file.")

        # Remove any chromosomes that aren't in the reference.
        input_variants_chroms = list(set(list(input_variants.CHROM)))
        for item in input_variants_chroms:
            if item not in reference_chromosomes and not printed_warning:
                print_and_log(f'Warning: ignoring all input vcf records for {item} '
                              f'because it is not found in the reference.', 'warning')
                print_and_log(f'\tIf this is unexpected, check that that {item} '
                              f'matches reference name exactly.', 'warning')
                printed_warning = True
                input_variants = input_variants[input_variants['CHROM'] != item]

        # Check the variants and classify as needed
        for chrom in reference_chromosomes:
            n_skipped = [0, 0, 0]
            if chrom in input_variants_chroms:
                for index, row in input_variants[input_variants['CHROM'] == chrom].iterrows():
                    span = (row['POS'], row['POS'] + len(row['REF']))
                    # -1 because going from VCF coords to array coords
                    r_seq = str(ref_index[chrom].seq[span[0] - 1:span[1] - 1])
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

                print_and_log(f'Found {len(input_variants)} valid variants for {chrom} in @include_vcf.', 'info')
                if any(n_skipped):
                    print_and_log(f'variants skipped: {sum(n_skipped)}', 'info')
                    print_and_log(f' - [{str(n_skipped[0])}] ref allele did not match reference', 'info')
                    print_and_log(f' - [{str(n_skipped[1])}] attempted to insert into N-region', 'info')
                    print_and_log(f' - [{str(n_skipped[2])}] alt allele contained non-ACGT characters', 'info')

    # parse input targeted regions, if present
    target_regions = parse_bed(options.target_bed, reference_chromosomes, begins_with_chr, False, options.debug)

    # parse discard bed similarly
    discard_regions = parse_bed(options.discard_bed, reference_chromosomes, begins_with_chr, False, options.debug)

    # parse input mutation rate rescaling regions, if present
    mutation_rate_regions, mutation_rate_values = parse_bed(options.mutation_bed, reference_chromosomes,
                                                            begins_with_chr, True, options.debug)

    """
    Initialize Output Files
    """
    if options.debug:
        logging.info("Beginning initialize output files.")

    # Prepare headers
    bam_header = None
    vcf_header = None
    if options.produce_bam:
        bam_header = ref_index
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
    if options.cancer:
        output_normal = out_prefix_parent_dir / f'{out_prefix_name}_normal'
        output_tumor = out_prefix_parent_dir / f'{out_prefix_name}_tumor'
        output_file_writer = OutputFileWriter(output_normal,
                                              paired=options.paired_ended,
                                              bam_header=bam_header,
                                              vcf_header=vcf_header,
                                              write_fastq=options.produce_fastq,
                                              write_fasta=options.produce_fasta,
                                              write_bam=options.produce_bam,
                                              write_vcf=options.produce_vcf
                                              )
        output_file_writer_cancer = OutputFileWriter(output_tumor,
                                                     paired=options.paired_ended,
                                                     bam_header=bam_header,
                                                     vcf_header=vcf_header,
                                                     write_fastq=options.produce_fastq,
                                                     write_fasta=options.produce_fasta,
                                                     write_bam=options.produce_bam,
                                                     write_vcf=options.produce_vcf
                                                     )
    else:
        outfile = out_prefix_parent_dir / out_prefix_name
        output_file_writer = OutputFileWriter(outfile,
                                              paired=options.paired_ended,
                                              bam_header=bam_header,
                                              vcf_header=vcf_header,
                                              write_fastq=options.produce_fastq,
                                              write_fasta=options.produce_fasta,
                                              write_bam=options.produce_bam,
                                              write_vcf=options.produce_vcf
                                              )

    print_and_log(f'Output files created.', 'info')

    """
    Begin Analysis
    """
    if options.debug:
        logging.info("Beginning Analysis")

    sequence = Seq("GATTACA")
    output_file_writer.write_fastq_record("bob", sequence, "BABBABA")

    cig = ['M'] * 101
    output_file_writer.write_bam_record(0, 'debug_neat-chrY-1', 2782028, cig, Seq('TTGGGAACACTCCTTAAACTTCAAAGAGAAGGCTCAAGGCCATTAAGAAAGACCTGTTGGGGTGCCATCCTTCCATCTGGCTTAGTCAGTTCTTAAAAAGA'), 'GCEGG6GDGG?FGG@\'%EGGGE@EGGE:FGEGFGGAGG?GGAGFFGGFGGGACFFGGDCGGC&/DG<EGGGG=-GEGGEG0DFF9F&D?G=:9DGFAG:6A', 99, 2782226, 70)

    output_file_writer.close_bam_file()

    # Step 1 (CAVA 469 - 475) - Find break points in the input file

    # Step 2 (CAVA 477 - 484) - Initialize simulation

    # Step 3 (CAVA 486 - 493) - Running simulation

    # Step 4 (CAVA 496 - 497) - Merging tmp files

    # Step 5 (CAVA 500 - 501) - Printing out summary and end time

    # We'll perform the analysis within a temp directory.
    with tempfile.TemporaryDirectory(prefix="sillyboy", dir=options.temp_dir) as temp_dir:
        pass

    print_and_log(f'NEAT completed in {datetime.datetime.now() - starttime}', 'info')
    print_and_log("Have a nice day!", 'info')
    print('-------------------------------------------------------------------\n')


if __name__ == '__main__':
    main()
