"""
The main class for parsing the input config file. This class reads in tho config yaml file using pyyaml,
Then updates any input fields. It also sets up the rng for the run. Because we wanted replicablity in this function,
we use a randomly generated seed, if no seed was input by the user. This is because numpy changed their
random number generator and you basically cannot retrieve the seed. See:
(https://stackoverflow.com/questions/32172054/how-can-i-retrieve-the-current-seed-of-numpys-random-number-generator)

So we pick a very large random number to use as a seed. The function will log this number in the log file,
and if you want to reproduce a run, you can input the seed in your config and you should get an exactly identical
output.

We want to introduce a quick-run option as well that will allow command line only inputs instead of a config
file. To do so, we have a few possible inputs the init function can accept. Inputting any of these will
trigger quick-run mode, setting most of the parameters to defaults.
"""

from types import SimpleNamespace
import numpy as np
import logging
import yaml
import importlib

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from typing import Final
from pathlib import Path

from ...common import validate_input_path


_LOG = logging.getLogger(__name__)


class Options(SimpleNamespace):
    """
    class representing the options
    """
    def __init__(self,
                 output_path: str | Path,
                 config_file: str | Path):
        SimpleNamespace.__init__(self)

        self.defs = {}
        self.config_file = config_file

        """
        Options flags for gen_reads. This metadata dict gives the type of variable (matching the python types)
        the default value ('.' means no default), and checks. There are four fields: option type (corresponding to
        a python type), default (or None, if no default), criteria 1 and criteria 2. Any items without values 
        should use None as a placeholder.

        Option Type: Current possible value are string, int, float, and boolean. These correspond to str, int, 
        float, and bool Python types. For files and directories, use string type.

        Criteria 1: There are two modes of checking: files and numbers. For files, criteria 1 should be set to
        'exists' to check file existence or None to skip a check, because that is not a user 
        input. For numbers, criteria 1 should be the lowest acceptable value (inclusive) for that variable.

        Criteria 2: For files, criteria 2 will not be checked, so set to None for consistency. For numbers, this
        should be the highest acceptable value (inclusive).
        (type, default, criteria1 (low/'exists'), criteria2 (high/None))
        """
        arbitrarily_large_number = 1e8
        self.defs['reference'] = ('string', None, 'exists', None)
        self.defs['partition_mode'] = ('string', 'chrom', None, None)
        self.defs['read_len'] = ('int', 101, 10, arbitrarily_large_number)
        self.defs['threads'] = ('int', 1, 1, arbitrarily_large_number)
        self.defs['coverage'] = ('float', 10.0, 1, arbitrarily_large_number)
        self.defs['error_model'] = ('string', None, 'exists', None)
        self.defs['avg_seq_error'] = ('float', None, 0, 0.3)
        self.defs['rescale_qualities'] = ('boolean', False, None, None)
        self.defs['ploidy'] = ('int', 2, 1, 100)
        self.defs['include_vcf'] = ('string', None, 'exists', None)
        self.defs['target_bed'] = ('string', None, 'exists', None)
        self.defs['discard_bed'] = ('string', None, 'exists', None)
        self.defs['off_target_scalar'] = ('float', 0.02, 0, 1)
        self.defs['discard_offtarget'] = ('boolean', False, None, None)
        self.defs['mutation_model'] = ('string', None, 'exists', None)
        self.defs['mutation_rate'] = ('float', None, 0, 0.3)
        self.defs['mutation_bed'] = ('string', None, 'exists', None)
        self.defs['n_handling'] = ('string', None, None, None)

        # Params for cancer (not implemented yet)
        self.defs['cancer'] = ('boolean', False, None, None)
        self.defs['cancer_model'] = ('string', None, 'exists', None)
        self.defs['cancer_purity'] = ('float', 0.8, 0.0, 1.0)

        self.defs['gc_model'] = ('string', None, 'exists', None)
        self.defs['paired_ended'] = ('boolean', False, None, None)
        self.defs['fragment_model'] = ('string', None, 'exists', None)
        self.defs['fragment_mean'] = ('float', None, 1, arbitrarily_large_number)
        self.defs['fragment_st_dev'] = ('float', None, 1, arbitrarily_large_number)
        self.defs['produce_bam'] = ('boolean', False, None, None)
        self.defs['produce_vcf'] = ('boolean', False, None, None)
        self.defs['produce_fasta'] = ('boolean', False, None, None)
        self.defs['output_config'] = ('boolean', False, None, None)
        self.defs['produce_fastq'] = ('boolean', True, None, None)
        self.defs['no_coverage_bias'] = ('boolean', False, None, None)
        self.defs['rng_seed'] = ('int', None, None, None)
        self.defs['min_mutations'] = ('int', 1, None, None)
        self.defs['fasta_per_ploid'] = ('boolean', False, None, None)

        # Cancer options (not yet implemented)
        self.cancer = False
        self.cancer_model = None
        self.cancer_purity = 0.8

        # Read the config file
        self.args = {}
        self.read()

        # Anything remaining is set to default:
        for key, (_, default, criteria1, criteria2) in self.defs.items():
            if key not in list(self.args.keys()):
                self.args[key] = default

        # Some options checking to clean up the args dict
        self.check_options()

        self.args['output'] = output_path

        self.__dict__.update(self.args)

        self.log_configuration()

    def set_value(self, key, value):
        if key in self.defs.keys():
            self.check_and_log_error(key, value, self.defs[key][2], self.defs[key][3])
        self.args[key] = value
        self.__dict__.update(self.args)

    @staticmethod
    def check_and_log_error(keyname, value_to_check, lowval, highval):
        if value_to_check is None:
            pass
        elif lowval != "exists" and highval:
            if not (lowval <= value_to_check <= highval):
                raise ValueError(f'@{keyname} must be between {lowval} and {highval}.')
        elif lowval == "exists" and value_to_check:
            validate_input_path(value_to_check, keyname)

    def read(self):
        config = yaml.load(open(self.config_file, 'r'), Loader=Loader)
        for key, value in config.items():
            if key in list(self.defs.keys()):
                type_of_var, default, criteria1, criteria2 = self.defs[key]
                # if it's already set to the default value, ignore.
                if value == default or value == ".":
                    continue

                # Now we check that the type is correct and it is in range, depending on the type defined for it
                # If it passes that it gets put into the args dictionary.
                if type_of_var == 'string':
                    temp = str(value)
                    self.check_and_log_error(key, temp, criteria1, criteria2)
                    self.args[key] = temp
                elif type_of_var == 'int':
                    try:
                        temp = int(value)
                    except Exception as ex:
                        _LOG.error(f'The value for {key} must be an integer. No decimals allowed.')
                        raise ex
                    self.check_and_log_error(key, temp, criteria1, criteria2)
                    self.args[key] = temp
                elif type_of_var == 'float':
                    try:
                        temp = float(value)
                    except ValueError:
                        _LOG.error(f'The value for {key} must be a float.')
                        raise
                    self.check_and_log_error(key, temp, criteria1, criteria2)
                    self.args[key] = temp
                elif type_of_var == 'boolean':
                    self.args[key] = value
                else:
                    raise KeyError(f'BUG: Undefined type in the Options dictionary: {type_of_var}.')

    def check_options(self):
        """
        Some sanity checks and corrections to the options.
        """

        # initialize random seed
        if not self.args['rng_seed']:
            """
            We want to allow for reproducibility in NEAT, so we'll need to generate a random number to set as the seed.
            I know this is statistically not as good, but Numpy made it more difficult to achieve a different way. 
            So we set the rng to a random state, then use that to pick a random integer, then use that as the seed for
            the simulation RNG.
            """
            seed_rng = np.random.default_rng()
            self.set_value('rng_seed', seed_rng.integers(2**52, 2**53, dtype=int))

        # Create the rng for this run
        self.set_value('rng', np.random.default_rng(self.args['rng_seed']))

        if not self.args['produce_bam'] and not self.args['produce_vcf'] \
                and not self.args['produce_fasta'] and not self.args['produce_fastq']:
            raise ValueError('No files would be produced, as all file types are set to false')

        # This next section just checks all the paired ended stuff
        flagged = False
        if self.args['paired_ended']:
            if self.args['fragment_model']:
                self.set_value('fragment_mean', None)
                self.set_value('fragment_st_dev', None)
            elif self.args['fragment_mean']:
                if not self.args['fragment_st_dev']:
                    flagged = True
            else:
                flagged = True
        if flagged:
            raise ValueError("For paired ended mode, you need to supply either a "
                             "@fragment_model or both @fragment_mean and @fragment_st_dev")

        self.set_value('n_handling', 'ignore')
        if self.args['paired_ended']:
            self.set_value('n_handling', 'random')

        # If discard_offtarget set to true and there is a targeted regions bed, set off_target_scalar to 0
        # If there is no targeted regions bed and discard_offtarget set to true, throw an error
        if self.args['discard_offtarget'] and self.args['target_bed']:
            self.set_value('off_target_scalar', 0.0)
        elif self.args['discard_offtarget'] and not self.args['target_bed']:
            _LOG.warning("@discard_offtarget set to true, but there is no target bed.")

    def log_configuration(self):
        """
        Combines the relevant parts of the input args and the options file to log a
        list of the configuration parameters. Useful for reproducibility.
        """
        _LOG.info(f'Run Configuration...')
        _LOG.info(f'Input fasta: {self.reference}')
        potential_filetypes = ['vcf', 'bam', 'fasta', 'fastq']
        extensions = []
        for suffix in potential_filetypes:
            key = f'produce_{suffix}'
            if self.args[key]:
                extensions.append(suffix)
        for item in extensions:
            _LOG.info(f'Output files: {self.output}.{item}')
        if self.threads == 1:
            _LOG.info(f"Single threading - 1 thread.")
        else:
            _LOG.warning(f'Multithreading coming soon!!')
            _LOG.info(f"Single threading - 1 thread.")
            # We'll work on multithreading later...
            # _LOG.info(f'Multithreading - {options.threads} threads')
        if self.paired_ended:
            _LOG.info(f'Running in paired-ended mode.')
            if self.fragment_model:
                _LOG.info(f'Using fragment length model: {self.fragment_model}')
            else:
                _LOG.info(f'Generating fragment model based on mean={self.fragment_mean}, '
                          f'st dev={self.fragment_st_dev}')
        else:
            _LOG.info(f'Running in single-ended mode.')
        _LOG.info(f'Using a read length of {self.read_len}')
        _LOG.info(f'Average coverage: {self.coverage}')
        if self.error_model:
            _LOG.info(f'Using error model: {self.error_model}')
        else:
            _LOG.info(f'Using default error model.')
        if self.avg_seq_error:
            _LOG.info(f'User defined average sequencing error rate: {self.avg_seq_error}.')
            if self.rescale_qualities:
                _LOG.info(f'Quality scores will be rescaled to match avg seq error rate.')
        _LOG.info(f'Ploidy value: {self.ploidy}')
        if self.include_vcf:
            _LOG.info(f'Vcf of variants to include: {self.include_vcf}')
        if self.target_bed:
            _LOG.info(f'BED of regions to target: {self.target_bed}')
            _LOG.info(f'Off-target coverage rate: {self.off_target_scalar}')
            _LOG.info(f'Discarding off target regions: {self.discard_offtarget}')
        if self.discard_bed:
            _LOG.info(f'BED of regions to discard: {self.discard_bed}')
        if self.mutation_model:
            _LOG.info(f'Using mutation model in file: {self.mutation_model}')
        if self.mutation_rate:
            _LOG.info(f'Rescaling average mutation rate to: {self.mutation_rate}')
        if self.mutation_bed:
            _LOG.info(f'BED of mutation rates of different regions: {self.mutation_bed}')
        if self.gc_model:
            _LOG.info(f'Using GC model: {self.gc_model}')
        if self.no_coverage_bias:
            _LOG.info(f'Ignoring GC bias from coverage calculations.')
        _LOG.info(f'RNG seed value for run: {self.rng_seed}')
