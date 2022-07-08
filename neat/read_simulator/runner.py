"""
Runner for generate_reads task
"""

import logging
import pickle
import gzip

from Bio import SeqIO
from pathlib import Path

from ..common import validate_input_path, validate_output_path
from .utils import Options, parse_input_vcf, parse_bed
from ..models import MutationModel, SequencingErrorModel, FragmentLengthModel, GcModel
from ..models.default_cancer_mutation_model import *

__all__ = ["read_simulator_runner"]

_LOG = logging.getLogger(__name__)


def initalize_all_models(options: Options):
    """
    Helper function that initializes models for use in the rest of the program.
    This includes loading the model and attaching the rng for this run
    to each model, so we can perform the various methods.

    :param options: the options for this run
    """

    if options.mutation_model:
        mut_model = pickle.load(gzip.open(options.mutation_model))
    else:
        mut_model = MutationModel()
    # Set the rng for the mutation model
    mut_model.rng = options.rng

    cancer_model = None
    if options.cancer and options.cancer_model:
        cancer_model = pickle.load(gzip.open(options.cancer_model))
        # Set the rng for the cancer mutation model
        cancer_model.rng = options.rng
    elif options.cancer:
        cancer_model = MutationModel(avg_mut_rate=default_cancer_avg_mut_rate,
                                     homozygous_freq=default_cancer_homozygous_freq,
                                     insertion_chance=default_cancer_insertion_chance,
                                     deletion_chance=default_cancer_deletion_chance,
                                     trinuc_trans_matrices=default_cancer_trinuc_trans_matrices,
                                     trinuc_trans_bias=default_cancer_trinuc_trans_bias,
                                     insertion_lengths=default_cancer_insertion_lengths,
                                     insertion_weights=default_cancer_insertion_weights,
                                     deletion_lengths=default_cancer_deletion_lengths,
                                     deletion_weights=default_cancer_deletion_weights,
                                     is_cancer=True,
                                     rng=options.rng)

    _LOG.debug("Mutation models loaded", 'debug')

    # We need sequencing errors to get the quality score attributes, even for the vcf
    if options.error_model:
        error_model = pickle.load(gzip.open(options.error_model))
    else:
        error_model = SequencingErrorModel()
    # Set the rng for the sequencing error model
    error_model.rng = options.rng

    _LOG.debug('Sequencing error model loaded', 'debug')

    # initialize gc_model
    if options.gc_model:
        gc_model = pickle.load(gzip.open(options.gc_model))
    else:
        gc_model = GcModel()
    # Set the rng for the GC bias model
    gc_model.rng = options.rng

    _LOG.debug('GC Bias model loaded', 'debug')

    if options.fragment_model:
        fraglen_model = pickle.load(gzip.open(options.fragment_model))
    else:
        fraglen_model = FragmentLengthModel()
    # Set the rng for the fragment length model
    fraglen_model.rng = options.rng

    _LOG.debug("Fragment length model loaded", 'debug')

    return mut_model, cancer_model, error_model, gc_model, fraglen_model


def read_simulator_runner(config: str, output: str):
    """
    Run the generate_reads function, which generates simulated mutations in a dataset and corresponding files.

    :param config: This is a configuration file. Keys start with @ symbol. Everything else is ignored.
    :param output: This is the prefix for the output.

    Raises
    ------
    FileNotFoundError
        If the given target or query file does not exist (or permissions
        prevent testing existence).
    RuntimeError
        If given target or query file exists but is empty.
    ValueError
        If neither prefix nor output path was provided.
    FileExistsError
        If the output file already exists.
    """
    _LOG.debug(f'config = {config}')
    _LOG.debug(f'output = {output}')

    validate_input_path(config)
    _LOG.info(f'Using configuration file {config}')
    validate_output_path(output, False)
    _LOG.info(f"Saving output files to {Path(output).parent}")

    options = Options(Path(output), config)

    """
    Model preparation

    Read input models or default models, as specified by user.
    """
    _LOG.info("Reading Models...")

    (
        mut_model,
        cancer_model,
        seq_error_model,
        gc_bias_model,
        fraglen_model
    ) = initalize_all_models(options)

    """
    Process Inputs
    """
    _LOG.info("Processing inputs.")
    _LOG.info(f'Reading {options.reference}.')

    reference_index = SeqIO.index(str(options.reference), 'fasta')
    _LOG.debug("Reference file indexed.")

    reference_contigs = list(reference_index.keys())
    options.set_value("reference_contigs", reference_contigs)

    input_variants = None
    if options.include_vcf:
        _LOG.info(f"Reading input VCF: {options.include_vcf}.")
        if options.cancer:
            (sample_names, input_variants) = parse_input_vcf(options.include_vcf,
                                                             options.ploidy,
                                                             mut_model.homozygous_freq,
                                                             reference_index,
                                                             options,
                                                             tumor_normal=True)

            tumor_ind = sample_names.index('tumor_sample')
            normal_ind = sample_names.index('normal_sample')
        else:
            (sample_names, input_variants) = parse_input_vcf(options.include_vcf,
                                                             options.ploidy,
                                                             mut_model.homozygous_freq,
                                                             reference_index,
                                                             options)

        _LOG.debug("Finished reading input vcf file")

    # parse input targeted regions, if present.
    if options.target_bed or options.discard_bed or options.mutation_bed:
        _LOG.info(f"Reading input bed files.")

    # Note that parse_bed will return None for any empty or missing files
    target_regions_dict = parse_bed(options.target_bed, options.reference_contigs,
                                    False)

    discard_regions_dict = parse_bed(options.discard_bed, options.reference_contigs,
                                     False)

    mutation_rate_dict = parse_bed(options.mutation_bed, options.reference_contigs,
                                   True)
