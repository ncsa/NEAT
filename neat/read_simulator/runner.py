"""
Runner for generate_reads task
"""

import logging
import pickle
import gzip

from Bio import SeqIO
from pathlib import Path

from .utils import Options, parse_input_vcf, parse_bed, OutputFileWriter, find_file_breaks, map_chromosome, \
    generate_variants
from ..common import validate_input_path, validate_output_path
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
    # Set user input mutation rate:
    if options.mutation_rate:
        mut_model.avg_mut_rate = options.mutation_rate

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

    _LOG.info(f'Using configuration file {config}')
    config = Path(config).resolve()
    validate_input_path(config)

    # prepare output
    _LOG.info(f"Saving output files to {Path(output).parent}")
    output = Path(output).resolve()

    if not output.parent.is_dir():
        _LOG.info('Creating output dir')
        output.parent.mkdir(parents=True, exist_ok=True)

    # Read options file
    options = Options(output, config)

    # Validate output
    validate_output_path(output, False, options.overwrite_output)

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

    # Note that parse_bed will return None for any empty or non-input files
    target_regions_dict = parse_bed(options.target_bed, options.reference_contigs,
                                    False)

    discard_regions_dict = parse_bed(options.discard_bed, options.reference_contigs,
                                     False)

    mutation_rate_dict = parse_bed(options.mutation_bed, options.reference_contigs,
                                   True)

    if any([target_regions_dict, discard_regions_dict, mutation_rate_dict]):
        _LOG.debug("Finished reading input beds.")

    # Prepare headers
    bam_header = None
    if options.produce_bam:
        # This is a dictionary that is the list of the contigs and the length of each.
        # This information will be needed later to create the bam header.
        bam_header = {key: len(reference_index[key]) for key in reference_index.keys()}

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf.
    # We'll also keep track here of what files we are producing.
    if options.cancer:
        output_normal = options.output.parent / f'{options.output.name}_normal'
        output_tumor = options.output.parent / f'{options.output.name}_tumor'
        output_file_writer = OutputFileWriter(options=options,
                                              bam_header=bam_header)
        output_file_writer_cancer = OutputFileWriter(options=options,
                                                     bam_header=bam_header)
    else:
        outfile = options.output.parent / options.output.name
        output_file_writer = OutputFileWriter(options=options,
                                              bam_header=bam_header)
        output_file_writer_cancer = None

    _LOG.debug(f'Output files ready for writing.')

    """
    Begin Analysis
    """
    _LOG.info("Beginning simulation.")

    breaks = find_file_breaks(options.threads, options.partition_mode, reference_index)

    _LOG.debug("Input reference partitioned for run")

    # these will be the features common to each contig, for multiprocessing
    common_features = {}

    vcf_files = []
    fasta_record = []
    fastq_files = []
    bam_files = []
    print_fasta_tell = False

    for contig in breaks:
        _LOG.info(f"Generating variants for {contig}")

        # Todo genericize breaks
        local_variants = None
        if input_variants:
            local_variants = {x[1]: input_variants[x] for x in input_variants if x[0] == contig}
        local_reference = reference_index[contig]

        _LOG.info(f'Creating trinucleotide map for {contig}...')
        local_trinuc_map = map_chromosome(local_reference, mut_model)

        # Since we're only running single threaded for now:
        threadidx = 1

        local_variant_file = options.temp_dir_path / f'{options.output.stem}_tmp_{contig}_{threadidx}.vcf'
        local_fasta_file = None
        if options.produce_fasta:
            local_fasta_file = options.temp_dir_path / f'{options.output.stem}_tmp_{contig}_{threadidx}.fasta'

        _LOG.debug(f'local vcf filename = {local_variant_file}')

        if threadidx == 1:
            # init_progress_info()
            pass

        generate_variants(reference=local_reference,
                          variant_file=local_variant_file,
                          trinucleotide_map=local_trinuc_map,
                          target_regions=target_regions_dict[contig],
                          discard_regions=discard_regions_dict[contig],
                          mutation_rate_regions=mutation_rate_dict[contig],
                          forced_variants=local_variants,
                          max_qual_score=max(seq_error_model.quality_scores),
                          mutation_model=mut_model,
                          options=options,
                          fasta_file=local_fasta_file)

        if options.produce_vcf:
            vcf_files.append(local_variant_file)

    if options.produce_vcf:
        _LOG.info("Outputting golden vcf.")
        output_file_writer.merge_temp_vcfs(vcf_files)

