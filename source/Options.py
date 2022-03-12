import os
import pathlib
from types import SimpleNamespace
import numpy as np
import random

from source.error_handling import premature_exit, print_and_log

# Constants
# NEAT Path is the path of the main NEAT executable, one level up from where this file resides
NEAT_PATH = pathlib.Path(__file__).resolve().parents[1]


class Options(SimpleNamespace):
    """
    class representing the options
    """

    def __init__(self, config_file=None):
        # For testing purposes, we'll allow a blank options object
        if not config_file:
            self.is_debug = True
        else:
            self.defs = {}
            self.config_file = config_file

            """
            Options flags for gen_reads. This metadata dict gives the type of variable (matching the python types)
            the default value ('.' means no default), and checks. There are four fields: option type (corresponding to
            a python type), default (or None, if no default), criteria 1 and criteria 2. Any items without values 
            should use None as a placeholder.

            Option Type: Current possible valeue are string, int, float, and boolean. These correspond to str, int, 
            float, and bool Python types. For files and directories, use string type.

            Criteria 1: There are two modes of checking: files and numbers. For files, criteria 1 should be set to
            'exists' to check file existence or None to skip a check (as with temp_dir, because that is not a user 
            input. For numbers, criteria 1 should be the lowest acceptable value (inclusive) for that variable.

            Criteria 2: For files, criteria 2 will not be checked, so set to None for consistency. For numbers, this
            should be the highest acceptable value (inclusive).
            (type, default, criteria1 (low/'exists'), criteria2 (high/None))
            """
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
            self.defs['temp_dir'] = ('string', pathlib.Path(os.getcwd()), None, None)
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
                                          f'(case-insensitive).', 'error')
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

        # initialize random seed
        if not self.args['rng_value']:
            self.args['rng_value'] = random.randint(2**52, 2**53)

        if not self.args['produce_bam'] and not self.args['produce_vcf'] \
                and not self.args['produce_fasta'] and not self.args['produce_fastq']:
            print_and_log('No files would be produced, as all file types are set to false', 'error')
            premature_exit(1)

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
            print_and_log("For paired ended mode, you need to supply either a "
                          "@fragment_model or both @fragment_mean and @fragment_st_dev", 'error')
            premature_exit(1)

        self.args['n_handling'] = "ignore"
        if self.args['paired_ended']:
            self.args['n_handling'] = "random"
