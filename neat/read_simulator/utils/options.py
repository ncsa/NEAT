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
import os
from copy import deepcopy
from typing import Any

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
                 reference: Path | None = None,
                 output_dir: Path | None = Path.cwd(),
                 fq1: Path | None = None,
                 fq2: Path | None = None,
                 vcf: Path | None = None,
                 bam: Path | None = None,
                 reads_pickle: Path | None = None,
                 output_prefix: str = "neat_sim",
                 overwrite_output: bool = False,
                 rng_seed: int | None = None,
                 read_len: int = 101,
                 threads: int = 1,
                 coverage: int = 10,
                 error_model: Path | None = None,
                 avg_seq_error: float | None = None,
                 rescale_qualities: bool = False,
                 ploidy: int = 2,
                 include_vcf: Path | None = None,
                 target_bed: Path | None = None,
                 discard_bed: Path | None = None,
                 mutation_model: Path | None = None,
                 mutation_rate: float | None = None,
                 mutation_bed: Path | None = None,
                 quality_offset: int = 33,
                 paired_ended: bool = False,
                 fragment_model: Path | None = None,
                 fragment_mean: float | None = None,
                 fragment_st_dev: float | None = None,
                 produce_bam: bool = False,
                 produce_vcf: bool = False,
                 produce_fastq: bool = True,
                 min_mutations: int = 0,
                 parallel_mode: str = "contig",
                 parallel_block_size: int = 500000,
                 cleanup_splits: bool = True,
                 splits_dir: Path | None = None,
                 reuse_splits: bool = False,
                 **kwargs: Any
                 ):
        """
        This is the base version of Options, which allows an empty reference,
        but we check this in the code to ensure we don't try to run without a
        reference.

        :param reference: The path to a reference file in fasta format
        :param output_dir: The path to the directory to write the output files
        :param fq1: The path with filename to the output fq1
        :param fq2: The path with filename to the output fq2
        :param vcf: The path with filename to the output vcf
        :param bam: the path with filename to the output bam
        :param output_prefix: The prefix to use for output files
        :param overwrite_output: If true, previous output will be overwritten, if
            filenames match
        :param rng_seed: The seed for reproducing these exact results
        :param read_len: The length of reads for the fastq output
        :param threads: The number of threads, or jobs, to run simultaneously.
        :param coverage: The target coverage depth
        :param error_model: A full path to a custom error model derived by model-sequencing-error
        :param avg_seq_error: A float to use for the target sequencing error. This must be between 0.0
            and 1.0
        :param rescale_qualities: To rescale quality scores to match the sequencing error
            model, set this to true.
        :param ploidy: Enter the number if the ploidy you wish to model does not equal 2
        :param include_vcf: Vcf of variants to include in the final output. Note that
            these must be of a simple form of ref and alt both being full sequences.
            i.e., in the case of a copy number 100 of variant ATG, Neat would need all
            100 ATGs written out to fully interpret
        :param target_bed: BED file of regions to target for mutations
        :param discard_bed: BED file of regions to discard
        :param mutation_model: Full path to model generated by gen-mut-model
        :param mutation_rate: Custom mutation rate for the dataset, overrides model mutation rate
        :param mutation_bed: BED indicating regions of interest and custom mutation rate for those
            regions
        :param quality_offset: Enter the number if the sequencing machine uses something other
            than PHRED+33 (enter number only, i.e., 33)
        :param paired_ended: Set to true to engage NEAT paired-ended mode
        :param fragment_model: Full path to model generated by model-fragment-lengths
        :param fragment_mean: To create a Normal distribution of fragment lengths, enter mean here
        :param fragment_st_dev: Same as above, but enter standard deviation for dataset here.
        :param produce_bam: True to create a bam file of alignments
        :param produce_vcf: True to produce a VCF of neat-added variants
        :param produce_fastq: False to turn off default fastq creation
        :param min_mutations: If you wish to have a minimunm number of mutations per block, enter it here
        :param parallel_mode: If you wish to use block size method, enter 'size' here
        :param parallel_block_size: If you use size method, specify any value but 500000 to change the block size
        :param cleanup_splits: Set to False in order to preserve splits after run
        :param reuse_splits: Attempts to reuse existing splits file
        """
        super().__init__(**kwargs)
        self.reference: Path = reference
        self.output_dir: Path = output_dir
        self.output_prefix: str = output_prefix
        self.overwrite_output = overwrite_output
        self.rng_seed = rng_seed
        self.read_len: int = read_len
        self.threads: int = threads
        self.coverage: int = coverage
        self.error_model: Path | None = error_model
        self.avg_seq_error: float | None = avg_seq_error
        self.rescale_qualities: bool = rescale_qualities
        self.ploidy: int = ploidy
        self.include_vcf: Path | None = include_vcf
        self.target_bed: Path | None = target_bed
        self.discard_bed: Path | None = discard_bed
        self.mutation_model: Path | None = mutation_model
        self.mutation_rate: float | None = mutation_rate
        self.mutation_bed: str | None = mutation_bed
        self.quality_offset: int = quality_offset

        self.paired_ended: bool = paired_ended
        self.fragment_model: str | None = fragment_model
        self.fragment_mean: float | None = fragment_mean
        self.fragment_st_dev: float | None = fragment_st_dev
        self.produce_bam: bool = produce_bam
        self.produce_vcf: bool = produce_vcf
        self.produce_fastq: bool = produce_fastq
        self.min_mutations: int = min_mutations

        # Parallel features
        self.parallel_mode: str = parallel_mode
        self.parallel_block_size: int = parallel_block_size
        self.cleanup_splits: bool = cleanup_splits
        self.splits_dir: Path | None = splits_dir
        self.reuse_splits: bool = reuse_splits

        # Actual output files
        self.fq1: Path | None = fq1
        self.fq2: Path | None = fq2
        self.vcf: Path | None = vcf
        self.bam: Path | None = bam
        # kind of a temporary holding pen for bam files
        self.reads_pickle: Path | None = reads_pickle

        # Set the rng for the run
        self.rng = self.set_random_seed()

        self.temporary_dir = TemporaryDirectory()
        self.temp_dir_path = Path(self.temporary_dir.name)

    @staticmethod
    def from_cli(output_dir: Path,
                 output_prefix: str,
                 config_file: Path,
                 ):
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
        # These are for checking
        defs = {
            'reference': (Path, None, 'exists', None),
            'read_len': (int, 151, 10, 1000000),
            'coverage': (int, 10, 1, 1000000),
            'error_model': (Path, None, 'exists', None),
            'avg_seq_error': (float, None, 0, 0.3),
            'rescale_qualities': (bool, False, None, None),
            'ploidy': (int, 2, 1, 100),
            'include_vcf': (Path, None, 'exists', None),
            'target_bed': (Path, None, 'exists', None),
            'discard_bed': (Path, None, 'exists', None),
            'mutation_model': (Path, None, 'exists', None),
            'mutation_rate': (float, None, 0, 0.3),
            'mutation_bed': (Path, None, 'exists', None),
            'quality_offset': (int, 33, 33, 64),
            'paired_ended': (bool, False, None, None),
            'fragment_model': (Path, None, 'exists', None),
            'fragment_mean': (float, None, 10, inf),
            'fragment_st_dev': (float, None, 1e-10, inf),
            'produce_bam': (bool, False, None, None),
            'produce_vcf': (bool, False, None, None),
            'produce_fastq': (bool, True, None, None),
            'no_coverage_bias': (bool, False, None, None),
            'rng_seed': (int, None, None, None),
            'min_mutations': (int, 0, None, None),
            'overwrite_output': (bool, False, None, None),
            'parallel_mode': (bool, False, None, None),
            'mode': (str, 'size', 'choice', ['size', 'contig']),
            'size': (int, 500000, None, None),
            'threads': (int, 1, 1, 1000),
            'cleanup_splits': (bool, True, None, None),
            'reuse_splits': (bool, False, None, None)
        }

        input_args = {}
        # Set up the dictionary of defaults
        for key, (_, default, _, _) in defs.items():
            if key not in input_args:
                input_args[key] = default

        base_options = Options(
            output_dir=output_dir,
            output_prefix=output_prefix
        )

        # Read the config file using the definitions for validation
        base_options.read_yaml(config_file, defs)

        # Merge validated config values over defaults
        final_args = dict(input_args)
        for key, val in defs.items():
            if not isinstance(val, tuple):
                final_args[key] = val

        # Update items to config or default values
        base_options.__dict__.update(final_args)
        base_options.set_random_seed()

        # Some options checking to clean up the args dict
        base_options.check_options()

        # Options not set by the config
        base_options.temporary_dir = TemporaryDirectory()
        base_options.temp_dir_path = Path(base_options.temporary_dir.name)

        base_options.log_configuration()
        return base_options

    @staticmethod
    def check_and_log_error(keyname, value_to_check, crit1, crit2):
        if value_to_check is None:
            pass
        elif crit1 == "exists" and value_to_check:
            validate_input_path(value_to_check)
        elif crit1 == "choice" and crit2:
            if value_to_check not in crit2:
                _LOG.error(f"Must choose one of {crit2}")
                sys.exit(1)
        elif isinstance(crit1, int) and isinstance(crit2, int):
            if not (crit1 <= value_to_check <= crit2):
                _LOG.error(f'`{keyname}` must be between {crit1} and {crit2} (input: {value_to_check}).')
                sys.exit(1)

    def read_yaml(self, config_yaml: Path, args: dict):
        """
        This sets up the option attributes. It's not perfect, because it sort of kills type hints.
        But I'm not sure how else to accomplish this.
        """
        config = yaml.load(open(config_yaml, 'r'), Loader=Loader)
        for key, value in config.items():
            if key in args:
                type_of_var, default, criteria1, criteria2 = args[key]
                # if it's already set to the default value, ignore.
                if value == default or value == ".":
                    continue
                # check for empty
                if value is None:
                    if key == "reference":
                        _LOG.error("Must entered a value for `reference` in config")
                        sys.exit(1)
                    else:
                        _LOG.debug(f"No value entered for `{key}`, skipping.")
                        continue

                # Now we check that the type is correct, and it is in range, depending on the type defined for it
                # If it passes that it gets put into the args dictionary.
                if type_of_var == Path:
                    if value != str(value):
                        _LOG.error(f"Incorrect type for value entered for {key}: {type_of_var} (found: {value})")
                        sys.exit(1)
                else:
                    if value != type_of_var(value):
                        _LOG.error(f"Incorrect type for value entered for {key}: {type_of_var} (found: {value})")
                        sys.exit(1)

                self.check_and_log_error(key, value, criteria1, criteria2)
                # Force reference to be a Path.
                if type_of_var == Path:
                    value = Path(value)
                args[key] = value

    def copy_with_changes(self,
                          reference: Path | None = None,
                          current_output_dir: Path | None = None,
                          fq1: Path | None = None,
                          fq2: Path | None = None,
                          reads_pickle: Path | None = None,
                          ):
        return_options = deepcopy(self)
        if reference is not None:
            return_options.reference = reference
        if current_output_dir is not None:
            return_options.current_output_dir = current_output_dir
        if fq1 is not None:
            return_options.fq1 = fq1
        if fq2 is not None:
            return_options.fq2 = fq2
        if reads_pickle is not None:
            return_options.reads_pickle = reads_pickle
            self.bam = None
        return return_options

    def set_random_seed(self) -> Generator:
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
        return np.random.default_rng(self.rng_seed)

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
                # Just turn these off so we don't try to load conflicting models
                self.fragment_mean = None
                self.fragment_st_dev = None

    def log_configuration(self):
        """
        Combines the relevant parts of the input args and the options file to log a
        list of the configuration parameters. Useful for reproducibility.
        """
        _LOG.info(f'Run Configuration...')
        _LOG.info(f'Input fasta: {self.reference}')
        _LOG.info(f'Outputting files to {self.output_dir}')
        files_to_produce = f'Producing the following files:\n'
        if self.produce_fastq:
            if self.paired_ended:
                fq1 = f'{str(self.output_dir)}/{self.output_prefix}_r1.fastq.gz'
                fq2 = f'{str(self.output_dir)}/{self.output_prefix}_r2.fastq.gz'
                validate_output_path(fq1, True, self.overwrite_output)
                validate_output_path(fq2, True, self.overwrite_output)
                files_to_produce += f'\t- {fq1}\n'
                files_to_produce += f'\t- {fq2}\n'
                self.fq1 = Path(fq1)
                self.fq2 = Path(fq2)
            else:
                fq1 = f'{str(self.output_dir)}/{self.output_prefix}.fastq.gz'
                validate_output_path(fq1, True, self.overwrite_output)
                files_to_produce += f'\t- {fq1}\n'
                self.fq1 = Path(fq1)
        if self.produce_bam:
            bam = f'{str(self.output_dir)}/{self.output_prefix}_golden.bam'
            validate_output_path(bam, True, self.overwrite_output)
            files_to_produce += f'\t- {bam}\n'
            self.bam = Path(bam)
        if self.produce_vcf:
            vcf = f'{str(self.output_dir)}/{self.output_prefix}_golden.vcf.gz'
            validate_output_path(vcf, True, self.overwrite_output)
            files_to_produce += f'\t- {vcf}\n'
            self.vcf = Path(vcf)

        _LOG.info(files_to_produce)

        if self.threads > 1:
            _LOG.info(f'Running read simulator in parallel mode.')
            _LOG.info(f'Multithreading - {self.threads} threads (or CPU Max)')
        else:
            _LOG.info(f"Single threading - 1 thread.")
            self.parallel_mode = 'contig'

        if self.parallel_mode == 'size':
            _LOG.info(f'Splitting reference into chunks.')
            _LOG.info(f'  - splitting input into size {self.size}')
        elif self.parallel_mode == 'contig':
            _LOG.info(f'Splitting input by contig.')
        if not self.cleanup_splits or self.reuse_splits:
            splits_dir = Path(f'{self.output_dir}/splits/')
            if splits_dir.is_dir():
                _LOG.info(f'Reusing existing splits {splits_dir}.')
            else:
                _LOG.warning(f'Reused splits set to True, but splits dir not found: {splits_dir}. Creating new splits')
            _LOG.info(f'Preserving splits for next run in directory {self.splits_dir}.')
        else:
            splits_dir = self.temp_dir_path / "splits"

        validate_output_path(splits_dir, False)
        self.splits_dir = splits_dir

        if self.produce_bam:
            try:
                import pysam
                _LOG.info(f"Using pysam: {pysam.__version__}")
            except ImportError:
                raise ImportError (
                    "Parallel NEAT requires pysam to be installed when produce_bam is set to true."
                )

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
