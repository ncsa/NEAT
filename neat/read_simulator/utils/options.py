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

import numpy as np
import logging
import yaml
import sys

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from types import SimpleNamespace
from tempfile import TemporaryDirectory
from pathlib import Path
from numpy.random import Generator
from math import inf

from ...common import validate_input_path, validate_output_path


_LOG = logging.getLogger(__name__)


class Options(SimpleNamespace):
    """
    class representing the options

    Note that for debugging/testing runs, we allow output_path and config_file to be empty.
    :param output_path: Path to the final output
    :param config_file: Path to the configuration file
    :param reference: Path to a reference for a test run
    :param rng_seed: A seed to use for a test run (for reproducible tests)
    """
    def __init__(self,
                 output_path: Path = None,
                 config_file: Path = None,
                 reference: Path = None,
                 rng_seed: int = None
                 ):

        SimpleNamespace.__init__(self)

        self.test_run = False
        if not output_path or not config_file:
            self.test_run = True

        """
        Options definitions for gen_reads. This metadata dict gives the type of variable (matching the python types)
        the default value ('.' means no default), and checks. There are four fields: option type (corresponding to
        a python type), default (or None, if no default), criteria 1 and criteria 2. Any items without values 
        should use None as a placeholder.

        Option Type: Current possible value are str, int, float, and bool.

        Criteria 1: There are two modes of checking: files and numbers. For files, criteria 1 should be set to
        'exists' to check file existence or None to skip a check, because that is not a user 
        input. For numbers, criteria 1 should be the lowest acceptable value (inclusive) for that variable.

        Criteria 2: For files, criteria 2 will not be checked, so set to None for consistency. For numbers, this
        should be the highest acceptable value (inclusive).
        (type, default, criteria1 (low/'exists'), criteria2 (high/None)
        """
        # TODO maybe redo as a dataclass?
        self.defs = {}
        self.config_file = config_file
        self.output = output_path

        self.defs['reference'] = (str, reference, 'exists', None)
        self.defs['read_len'] = (int, 151, 10, inf)
        self.defs['threads'] = (int, 1, 1, inf)
        self.defs['coverage'] = (int, 10, 1, inf)
        self.defs['error_model'] = (str, None, 'exists', None)
        self.defs['avg_seq_error'] = (float, None, 0, 0.3)
        self.defs['rescale_qualities'] = (bool, False, None, None)
        self.defs['ploidy'] = (int, 2, 1, 100)
        self.defs['include_vcf'] = (str, None, 'exists', None)
        self.defs['target_bed'] = (str, None, 'exists', None)
        self.defs['discard_bed'] = (str, None, 'exists', None)
        self.defs['mutation_model'] = (str, None, 'exists', None)
        self.defs['mutation_rate'] = (float, None, 0, 0.3)
        self.defs['mutation_bed'] = (str, None, 'exists', None)
        self.defs['quality_offset'] = (int, 33, 33, 64)

        # Params for cancer (not implemented yet)
        self.defs['cancer'] = (bool, False, None, None)
        self.defs['cancer_model'] = (str, None, 'exists', None)
        self.defs['cancer_purity'] = (float, 0.8, 0.0, 1.0)

        self.defs['paired_ended'] = (bool, False, None, None)
        self.defs['fragment_model'] = (str, None, 'exists', None)
        self.defs['fragment_mean'] = (float, None, 10, inf)
        self.defs['fragment_st_dev'] = (float, None, 1e-10, inf)
        self.defs['produce_bam'] = (bool, False, None, None)
        self.defs['produce_vcf'] = (bool, False, None, None)
        self.defs['produce_fastq'] = (bool, True, None, None)
        self.defs['no_coverage_bias'] = (bool, False, None, None)

        # These are primarily debug options
        self.defs['rng_seed'] = (int, None, None, None)
        self.defs['min_mutations'] = (int, 0, None, None)
        self.defs['overwrite_output'] = (bool, False, None, None)

        # Create base variables, for update by the config
        self.reference: str | Path = reference
        self.read_len: int = 151
        self.threads: int = 1
        self.coverage: int = 10
        self.error_model: str | None = None
        self.avg_seq_error: float | None = None
        self.rescale_qualities: bool = False
        self.ploidy: int = 2
        self.include_vcf: str | None = None
        self.target_bed: str | None = None
        self.discard_bed: str | None = None
        self.mutation_model: str | None = None
        self.mutation_rate: float | None = None
        self.mutation_bed: str | None = None
        self.quality_offset: int = 33

        self.paired_ended: bool = False
        self.fragment_model: str | None = None
        self.fragment_mean: float | None = None
        self.fragment_st_dev: float | None = None
        self.produce_bam: bool = False
        self.produce_vcf: bool = False
        self.produce_fastq: bool = True

        # These are primarily debug options.
        self.min_mutations: int = 0
        self.overwrite_output: bool = False
        self.rng_seed: int | None = None
        self.rng: Generator | None = None

        # Cancer options (not yet implemented)
        self.cancer: str | Path
        self.cancer_model: bool
        self.cancer_purity: float

        self.args = {}
        # Set up the dictionary
        for key, (_, default, _, _) in self.defs.items():
            if key not in self.args:
                self.args[key] = default

        # Read the config file
        if not self.test_run:
            self.read()

        # Update items to config or default values
        self.__dict__.update(self.args)

        if self.test_run:
            self.rng_seed = rng_seed
        # Set the rng for the run
        self.set_random_seed()

        # Some options checking to clean up the args dict
        if not self.test_run:
            self.check_options()

        # Options not set by the config
        self.temporary_dir = TemporaryDirectory()
        self.temp_dir_path = Path(self.temporary_dir.name)
        # Set later after reference processing
        self.reference_contigs = None

        if not self.test_run:
            self.log_configuration()

    @staticmethod
    def check_and_log_error(keyname, value_to_check, lowval, highval):
        if value_to_check is None:
            pass
        elif lowval != "exists" and highval:
            if not (lowval <= value_to_check <= highval):
                _LOG.error(f'`{keyname}` must be between {lowval} and {highval} (input: {value_to_check}).')
                sys.exit(1)

        elif lowval == "exists" and value_to_check:
            validate_input_path(value_to_check)

    def read(self):
        """
        This sets up the option attributes. It's not perfect, because it sort of kills type hints.
        But I'm not sure how else to accomplish this.
        """
        # Skip trying to read the config for a test run
        config = yaml.load(open(self.config_file, 'r'), Loader=Loader)
        for key, value in config.items():
            if key in self.defs:
                type_of_var, default, criteria1, criteria2 = self.defs[key]
                # if it's already set to the default value, ignore.
                if value == default or value == ".":
                    continue
                # check for empty
                if value is None:
                    if key == "reference":
                        _LOG.error("Must entered a value for `reference` in config")
                        sys.exit(1)
                    else:
                        _LOG.warning(f"No value entered for `{key}`, skipping.")
                        continue

                # Now we check that the type is correct, and it is in range, depending on the type defined for it
                # If it passes that it gets put into the args dictionary.
                if value != type_of_var(value):
                    _LOG.error(f"Incorrect type for value entered for {key}: {type_of_var} (found: {value})")
                    sys.exit(1)

                self.check_and_log_error(key, value, criteria1, criteria2)
                self.args[key] = value

    def set_random_seed(self):
        """
        Sets up random number generator, which will be used for the run.
        """
        # 0 is a valid entry, but it fails this check with just "not self.rng_seed", so we have to be explicit
        if not self.rng_seed and not self.rng_seed == 0:
            """
            We want to allow for reproducibility in NEAT, so we'll need to generate a random number to set as the seed.
            I know this is statistically not as good, but Numpy made it more difficult to achieve a different way. 
            So we set the rng to a random state, then use that to pick a random integer, then use that as the seed for
            the simulation RNG.
            """
            seed_rng = np.random.default_rng()
            # I pulled this idea from the internet re 2^52 - 2^53. Supposed to more reliable for randomness
            self.rng_seed = seed_rng.integers(2 ** 52, 2 ** 53, dtype=int)

        # Create the rng for this run
        self.rng = np.random.default_rng(self.rng_seed)

    def check_options(self):
        """
        Some sanity checks and corrections to the options.
        """
        if not (self.produce_bam or self.produce_vcf or self.produce_fastq):
            _LOG.error('No files would be produced, as all file types are set to false')
            sys.exit(1)

        # This next section just checks all the paired ended stuff
        if self.paired_ended:
            if self.fragment_model:
                self.fragment_mean = None
                self.fragment_st_dev = None

    def log_configuration(self):
        """
        Combines the relevant parts of the input args and the options file to log a
        list of the configuration parameters. Useful for reproducibility.
        """
        _LOG.info(f'Run Configuration...')
        _LOG.info(f'Input fasta: {self.reference}')
        files_to_produce = f'Producing the following files:\n'
        if self.produce_fastq:
            if self.paired_ended:
                fq1 = f'{self.output}_r1.fastq.gz'
                fq2 = f'{self.output}_r2.fastq.gz'
                validate_output_path(fq1, True, self.overwrite_output)
                validate_output_path(fq2, True, self.overwrite_output)
                files_to_produce += f'\t- {fq1}\n'
                files_to_produce += f'\t- {fq2}\n'
            else:
                fq1 = f'{self.output}.fastq.gz'
                validate_output_path(fq1, True, self.overwrite_output)
                files_to_produce += f'\t- {fq1}\n'
        if self.produce_bam:
            bam = f'{self.output}_golden.bam'
            validate_output_path(bam, True, self.overwrite_output)
            files_to_produce += f'\t- {bam}\n'
        if self.produce_vcf:
            vcf = f'{self.output}_golden.vcf.gz'
            validate_output_path(vcf, True, self.overwrite_output)
            files_to_produce += f'\t- {vcf}\n'

        _LOG.info(files_to_produce)

        if self.threads == 1:
            _LOG.info(f"Single threading - 1 thread.")
        else:
            _LOG.warning(f'Multithreading coming soon!!')
            _LOG.info(f"Single threading - 1 thread.")
            # We'll work on multithreading later...
            # _LOG.info(f'Multithreading - {options.threads} threads')
        _LOG.info(f'Using a read length of {self.read_len}')
        if self.fragment_mean:
            if self.fragment_mean < self.read_len:
                _LOG.error(f"`fragment_mean` (input: {self.fragment_mean}) "
                           f"must be at least as long as `read_len` (input or default: {self.read_len}).")
                sys.exit(1)
            if self.fragment_st_dev:
                _LOG.info(f'Generating fragments based on mean={self.fragment_mean}, '
                          f'stand. dev={self.fragment_st_dev}')
            else:
                _LOG.error("Provided fragment mean, but no fragment standard deviation!")
                sys.exit(1)
        if self.paired_ended:
            _LOG.info(f'Running in paired-ended mode.')
            if self.fragment_model:
                _LOG.info(f'Using fragment length model: {self.fragment_model}')
            elif self.fragment_mean:
                pass  # Already addressed this above
            else:
                _LOG.error("Paired ended mode requires either a fragment model or a mean/standard deviation.")
                sys.exit(1)
        else:
            _LOG.info(f'Running in single-ended mode.')
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
        if self.discard_bed:
            _LOG.info(f'BED of regions to discard: {self.discard_bed}')
        if self.mutation_model:
            _LOG.info(f'Using mutation model in file: {self.mutation_model}')
        if self.mutation_rate:
            _LOG.info(f'Custom average mutation rate for the run: {self.mutation_rate}')
        if self.mutation_bed:
            _LOG.info(f'BED of mutation rates of different regions: {self.mutation_bed}')
        _LOG.info(f'RNG seed value for run: {self.rng_seed}')
