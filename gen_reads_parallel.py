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
import gzip
import logging
import multiprocessing
from Bio import SeqIO
from types import SimpleNamespace
from source.ref_func import index_ref, read_ref

from gen_reads import main as processing_main


# Constants
VERSION = 3.0
NEAT_PATH = pathlib.Path(__file__).resolve().parent


# useful functions for the rest
def print_start_info() -> datetime.datetime:
    """
    Prints a NEAT start message with version number and logs the start time

    :return: logs the time this function was called as start time.
    """
    starttime = datetime.datetime.now()
    print('\n-----------------------------------------------------------')
    print(f"NEAT multithreaded version, v{VERSION}, is running.")
    print(f'Started: {str(starttime)}.')
    logging.info(f'NEAT multithreaded version, v{VERSION} started.')
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
                print(f'\nERROR: @{keyname} must be between {lowval} and {highval}.')
                print(f'\nNothing written, quitting NEAT.')
                print('-----------------------------------------------------------------------')
                logging.error(f'ERROR: @{keyname} must be between {lowval} and {highval}.')
                sys.exit(1)
        elif lowval == "exists":
            if not pathlib.Path(value_to_check).is_file():
                print(f'\nERROR: the file given to @{keyname} does not exist')
                print(f'\nNothing written, quitting NEAT.')
                print('-----------------------------------------------------------------------')
                logging.error(f'ERROR: the file given to @{keyname} does not exist')
                sys.exit(1)
        elif not lowval and not highval:
            # This indicates a boolean and we have nothing to check
            pass
        else:
            print(f'\nUndeclared criteria {lowval, highval} in Options definitions.')
            logging.error(f'Undeclared criteria {lowval, highval} in Options definitions.')
            sys.exit(1)

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
                            print(f'\nThe value for {key} must be an integer. No decimals allowed.')
                            logging.error(f'The value for {key} must be an integer. No decimals allowed.')
                            sys.exit(1)
                        self.check_and_log_error(key, temp, criteria1, criteria2)
                        self.args[key] = temp
                    elif type_of_var == 'float':
                        try:
                            temp = float(line_split[1])
                        except ValueError:
                            print(f'\nThe value for {key} must be a float.')
                            logging.error(f'The value for {key} must be a float.')
                            sys.exit(1)
                        self.check_and_log_error(key, temp, criteria1, criteria2)
                        self.args[key] = temp
                    elif type_of_var == 'boolean':
                        try:
                            if line_split[1].lower() == 'true':
                                self.args[key] = True
                            elif line_split.lower() == 'false':
                                self.args[key] = False
                            else:
                                raise ValueError
                        except ValueError:
                            print(f'\nBoolean key @{key} '
                                  f'requires a value of "true" or "false" (case insensitive).')
                            logging.error(f'Boolean key @{key} '
                                          f'requires  a value of "true" or "false" (case insensitive).')
                            sys.exit(1)
                    else:
                        print(f'\nBUG: Undefined type in the Options dictionary: {type_of_var}.')
                        logging.error(f'BUG: Undefined type in the Options dictionary: {type_of_var}.')
                        sys.exit(1)
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
        if self.args['threads'] == 1:
            print("\nRunning in single threaded mode")
            logging.info("Running in single threaded mode")
        if self.args['produce_fasta']:
            print("\nFASTA mode active.")
            print("NOTE: At the moment, NEAT can produce a FASTA or FASTQ files, not both.")
            logging.info("FASTA mode active.")
            logging.info("NOTE: At the moment, NEAT can produce a FASTA OR FASTQ files, not both.")
            # Turn off fastq and paired-ended mode for now
            self.args['produce_fastq'] = False
            self.args['paired_ended'] = False
            self.args['fragment_model'] = None
            self.args['fragment_mean'] = None
            self.args['fragment_st_dev'] = None
        if not self.args['produce_bam'] and not self.args['produce_vcf'] \
                and not self.args['produce_fasta'] and not self.args['produce_fastq']:
            print('\nERROR: No files would be produced, all file types turned off')
            logging.error('ERROR: No files would be produced, all file types turned off')
            sys.exit(1)
        if not self.args['produce_fastq']:
            print("\nBypassing FASTQ generation.")
            logging.info("Bypassing FASTQ generation.")
        if self.args['produce_vcf'] and (not self.args['produce_fastq'] and not self.args['produce_bam']
                                         and not self.args['produce_fasta']):
            print('Only producing VCF output.')
            logging.info('Only producing VCF output')
        if self.args['produce_bam'] and (not self.args['produce_fastq'] and not self.args['produce_vcf']
                                         and not self.args['produce_fasta']):
            print('Only producing BAM output.')
            logging.info('Only producing BAM output')

        # This next section just checks all the paired ended stuff
        flagged = False
        if self.args['paired_ended']:
            print("\nPaired-ended mode")
            logging.info('Paired-ended mode')
            if self.args['fragment_model']:
                print(f"\nUsing fragment length model {self.args['fragment_model']} to produce paired ended reads")
                logging.info(f"Using fragment length model {self.args['fragment_model']} to produce paired ended reads")
                self.args['fragment_mean'] = None
                self.args['fragment_st_dev'] = None
            elif self.args['fragment_mean']:
                if self.args['fragment_st_dev']:
                    print(f"\nUsing fragment length mean = {self.args['fragment_mean']}, "
                          f"standard deviation = {self.args['fragment_st_dev']} to produce paired ended reads.")
                    logging.info(f"Using fragment length mean = {self.args['fragment_mean']}, "
                                 f"standard deviation = {self.args['fragment_st_dev']} to produce paired ended reads.")
                else:
                    flagged = True
            else:
                flagged = True
        if flagged:
            print("\nERROR: For paired ended mode, you need either a "
                  "@fragment_model or both @fragment_mean and @fragment_st_dev")
            logging.error("ERROR: For paired ended mode, you need either a "
                          "@fragment_model or both @fragment_mean and @fragment_st_dev")
            sys.exit(1)


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

    I'm thinking we don't really need to parse other args for this, but just pass them along, but I could be wrong.
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

    starttime = print_start_info()

    log_name = f'{args.output}.log'

    logging.basicConfig(filename=log_name, filemode='w',
                        format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.DEBUG)

    if not args.conf.is_file():
        print(f'\nError: configuration file ({args.conf}) cannot be found.')
        logging.error(f'Error: configuration file ({args.conf}) cannot be found.')
        sys.exit(1)

    # Reads in the user-entered options from the config file and performs
    # some basic checks and corrections.
    user_options = Options(args.conf)

    ref_index = SeqIO.index(str(user_options.reference), 'fasta')
    reference_chromosomes = list(ref_index.keys())

    print(reference_chromosomes)
    print(f'It took {datetime.datetime.now() - starttime}')
    print("Have a nice day!")


if __name__ == '__main__':
    main()
