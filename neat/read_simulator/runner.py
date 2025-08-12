"""
Runner for generate_reads task
"""
import time
import logging
import pickle
import gzip

from Bio import SeqIO
from pathlib import Path

from .utils import Options, parse_input_vcf, parse_beds, OutputFileWriter, \
    generate_variants, generate_reads
from ..common import validate_input_path, validate_output_path
from ..models import MutationModel, SequencingErrorModel, FragmentLengthModel, TraditionalQualityModel
from ..models.default_cancer_mutation_model import *
from ..variants import ContigVariants

__all__ = ["read_simulator_runner"]

_LOG = logging.getLogger(__name__)


def initialize_all_models(options: Options):
    """
    Helper function that initializes models for use in the rest of the program.
    This includes loading the model and attaching the rng for this run
    to each model, so we can perform the various methods.

    :param options: the options for this run
    """

    # Load mutation model or instantiate default
    if options.mutation_model:
        mut_model = pickle.load(gzip.open(options.mutation_model))
    else:
        mut_model = MutationModel()

    # Set random number generator for the mutations:
    mut_model.rng = options.rng
    # Set custom mutation rate for the run, or set the option to the input rate so we can use it later
    if options.mutation_rate is not None:
        mut_model.avg_mut_rate = options.mutation_rate

    cancer_model = None
    if options.cancer and options.cancer_model:
        # cancer_model = pickle.load(gzip.open(options.cancer_model))
        # Set the rng for the cancer mutation model
        cancer_model.rng = options.rng
    elif options.cancer:
        # Note all parameters not entered here use the mutation madel defaults
        cancer_model = MutationModel(
            avg_mut_rate=default_cancer_avg_mut_rate,
            homozygous_freq=default_cancer_homozygous_freq,
            variant_probs=default_cancer_variant_probs,
            insert_len_model=default_cancer_insert_len_model,
            is_cancer=True
        )

    _LOG.debug("Mutation models loaded")

    # We need sequencing errors to get the quality score attributes, even for the vcf
    if options.error_model:
        error_models = pickle.load(gzip.open(options.error_model))
        error_model_1 = error_models["error_model1"]
        quality_score_model_1 = error_models["qual_score_model1"]
        if options.paired_ended:
            if error_models["error_model2"]:
                error_model_2 = error_models["error_model2"]
                quality_score_model_2 = error_models["qual_score_model2"]
            else:
                _LOG.warning('Paired ended mode declared, but input sequencing error model is single ended,'
                             'duplicating model for both ends')
                error_model_2 = error_models["error_model1"]
                quality_score_model_2 = error_models["qual_score_model1"]
        else:
            # ignore second model if we're in single-ended mode
            error_model_2 = None
            quality_score_model_2 = None
    else:
        # Use all the default values
        error_model_1 = SequencingErrorModel()
        quality_score_model_1 = TraditionalQualityModel()
        if options.paired_ended:
            error_model_2 = SequencingErrorModel()
            quality_score_model_2 = TraditionalQualityModel()
        else:
            error_model_2 = None
            quality_score_model_2 = None

    _LOG.debug('Sequencing error and quality score models loaded')

    if options.fragment_model:
        fraglen_model = pickle.load(gzip.open(options.fragment_model))
        fraglen_model.rng = options.rng
    elif options.fragment_mean:
        fraglen_model = FragmentLengthModel(options.fragment_mean, options.fragment_st_dev)
    else:
        # For single ended, fragment length will be based on read length
        fragment_mean = options.read_len * 2.0
        fragment_st_dev = fragment_mean * 0.2
        fraglen_model = FragmentLengthModel(fragment_mean, fragment_st_dev)

    _LOG.debug("Fragment length model loaded")

    return \
        mut_model, \
        cancer_model, \
        error_model_1, \
        error_model_2, \
        quality_score_model_1, \
        quality_score_model_2, \
        fraglen_model


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
        seq_error_model_1,
        seq_error_model_2,
        qual_score_model_1,
        qual_score_model_2,
        fraglen_model
    ) = initialize_all_models(options)

    """
    Process Inputs
    """
    _LOG.info(f'Reading {options.reference}.')

    # TODO check into SeqIO.index_db()
    reference_index = SeqIO.index(str(options.reference), "fasta")
    reference_keys_with_lens = {key: len(value) for key, value in reference_index.items()}
    _LOG.debug("Reference file indexed.")

    if _LOG.getEffectiveLevel() < 20:
        count = 0
        for contig in reference_keys_with_lens:
            count += reference_keys_with_lens[contig]
        _LOG.debug(f"Length of reference: {count / 1_000_000:.2f} Mb")

    input_variants_dict = {x: ContigVariants() for x in reference_keys_with_lens}
    if options.include_vcf:
        _LOG.info(f"Reading input VCF: {options.include_vcf}.")
        if options.cancer:
            # TODO Check if we need full ref index or just keys and lens
            sample_names = parse_input_vcf(
                input_variants_dict,
                options.include_vcf,
                options.ploidy,
                mut_model.homozygous_freq,
                reference_index,
                options,
                tumor_normal=True
            )

            tumor_ind = sample_names['tumor_sample']
            normal_ind = sample_names['normal_sample']
        else:
            # TODO Check if we need full ref index or just keys and lens
            sample_names = parse_input_vcf(
                input_variants_dict,
                options.include_vcf,
                options.ploidy,
                mut_model.homozygous_freq,
                reference_index,
                options
            )

        _LOG.debug("Finished reading input vcf file")

    # Note that parse_beds will return None for any empty or non-input files
    bed_files = (options.target_bed, options.discard_bed, options.mutation_bed)

    # parse input targeted regions, if present.
    if any(bed_files):
        _LOG.info(f"Reading input bed files.")

    (
        target_regions_dict,
        discard_regions_dict,
        mutation_rate_dict
    ) = parse_beds(options, reference_keys_with_lens, mut_model.avg_mut_rate)

    if any(bed_files):
        _LOG.debug("Finished reading input beds.")

    # Prepare headers
    bam_header = None
    if options.produce_bam:
        # This is a dictionary that is the list of the contigs and the length of each.
        # This information will be needed later to create the bam header.
        bam_header = reference_keys_with_lens

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

    breaks = find_file_breaks(reference_keys_with_lens)

    _LOG.debug("Input reference partitioned for run")

    # these will be the features common to each contig, for multiprocessing
    common_features = {}

    local_variant_files = {}
    fastq_files = []

    sam_reads_files = []

    for contig in breaks:
        local_variant_files[contig] = None

        _LOG.info(f"Generating variants for {contig}")

        input_variants = input_variants_dict[contig]
        # TODO: add the ability to pick up input variants here from previous loop

        local_reference = reference_index[contig]

        # Since we're only running single threaded for now:
        threadidx = 1

        local_bam_pickle_file = None
        if options.produce_bam:
            local_bam_pickle_file = options.temp_dir_path / f'{options.output.stem}_tmp_{contig}_{threadidx}.p.gz'

        if threadidx == 1:
            # init_progress_info()
            pass

        if options.paired_ended:
            max_qual_score = max(max(qual_score_model_1.quality_scores), max(qual_score_model_2.quality_scores))
        else:
            max_qual_score = max(qual_score_model_1.quality_scores)

        local_variants = generate_variants(
            reference=local_reference,
            mutation_rate_regions=mutation_rate_dict[contig],
            existing_variants=input_variants,
            mutation_model=mut_model,
            max_qual_score=max_qual_score,
            options=options
        )

        # This function saves the local variant data a dictionary. We may need to write this to file.
        local_variant_files[contig] = local_variants

        if options.produce_fastq or options.produce_bam:
            read1_fastq_paired, read1_fastq_single, read2_fastq_paired, read2_fastq_single = generate_reads(
                local_reference,
                local_bam_pickle_file,
                seq_error_model_1,
                seq_error_model_2,
                qual_score_model_1,
                qual_score_model_2,
                fraglen_model,
                local_variants,
                options.temp_dir_path,
                target_regions_dict[contig],
                discard_regions_dict[contig],
                options,
                contig,
                )

            contig_temp_fastqs = ((read1_fastq_paired, read2_fastq_paired), (read1_fastq_single, read2_fastq_single))
            fastq_files.append(contig_temp_fastqs)
            if options.produce_bam:
                sam_reads_files.append(local_bam_pickle_file)

    if options.produce_vcf:
        _LOG.info(f"Outputting golden vcf: {str(output_file_writer.vcf_fn)}")
        output_file_writer.write_final_vcf(local_variant_files, reference_index)

    if options.produce_fastq:
        if options.paired_ended:
            _LOG.info(f"Outputting fastq files: "
                      f"{', '.join([str(x) for x in output_file_writer.fastq_fns]).strip(', ')}")
        else:
            _LOG.info(f"Outputting fastq file: {output_file_writer.fastq_fns[0]}")
        output_file_writer.merge_temp_fastqs(fastq_files, options.rng)

    if options.produce_bam:
        _LOG.info(f"Outputting golden bam file: {str(output_file_writer.bam_fn)}")
        contig_list = list(reference_keys_with_lens)
        contigs_by_index = {contig_list[n]: n for n in range(len(contig_list))}
        output_file_writer.output_bam_file(sam_reads_files, contigs_by_index, options.read_len)


# Initial implementation of parallelization
def find_file_breaks(reference_keys_with_lens: dict) -> dict:
    """
    Returns a dictionary with the chromosomes as keys, which is the start of building the chromosome map

    :param reference_keys_with_lens: a dictionary with chromosome keys and sequence values
    :return: a dictionary containing the chromosomes as keys and either "all" for values, or a list of indices
    """
    partitions = {}
    for contig in reference_keys_with_lens:
        partitions[contig] = [(0, reference_keys_with_lens[contig])]

    return partitions

# TO DO: Add function to read in parallelization from config
