"""
Runner for generate_reads task
"""
import copy
import logging
import pickle
import gzip

from Bio import SeqIO
from pathlib import Path

from .utils import Options, parse_input_vcf, parse_beds, fill_out_bed_dict, OutputFileWriter, find_file_breaks, \
    map_chromosome, generate_variants, write_local_file, generate_reads
from ..common import validate_input_path, validate_output_path
from ..models import MutationModel, SequencingErrorModel, FragmentLengthModel, GcModel
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
    if options.mutation_rate:
        mut_model.avg_mut_rate = options.mutation_rate

    cancer_model = None
    if options.cancer and options.cancer_model:
        cancer_model = pickle.load(gzip.open(options.cancer_model))
        # Set the rng for the cancer mutation model
        cancer_model.rng = options.rng
    elif options.cancer:
        # Note all parameters not entered here use the mutation madel defaults
        cancer_model = MutationModel(avg_mut_rate=default_cancer_avg_mut_rate,
                                     homozygous_freq=default_cancer_homozygous_freq,
                                     variant_probs=default_cancer_variant_probs,
                                     insert_len_model=default_cancer_insert_len_model,
                                     is_cancer=True,
                                     rng=options.rng)

    _LOG.debug("Mutation models loaded")

    # We need sequencing errors to get the quality score attributes, even for the vcf
    if options.error_model:
        error_models = pickle.load(gzip.open(options.error_model))
        error_model_1 = error_models[0]
        if options.paired_ended:
            if error_models[1]:
                error_model_2 = error_models[1]
            else:
                _LOG.warning('Paired ended mode declared, but input sequencing error model is single ended,'
                             'duplicating model for both ends')
                error_model_2 = error_models[0]
        else:
            # ignore second model if we're in single-ended mode
            error_model_2 = None
    else:
        # Use all the default values
        error_model_1 = SequencingErrorModel()
        if options.paired_ended:
            error_model_2 = SequencingErrorModel()
        else:
            error_model_2 = None
    # Set the rng for the sequencing error model
    error_model_1.rng = options.rng
    if error_model_2:
        error_model_2.rng = options.rng

    _LOG.debug('Sequencing error models loaded')

    # initialize gc_model
    if options.gc_model:
        gc_model = pickle.load(gzip.open(options.gc_model))
    else:
        # Creates a default model
        gc_model = GcModel()
    # Set the coverage target for the GC bias model
    gc_model.coverage = options.coverage

    _LOG.debug('GC Bias model loaded')

    fraglen_model = None
    if options.paired_ended:
        if options.fragment_model:
            fraglen_model = pickle.load(gzip.open(options.fragment_model))
        else:
            fraglen_model = FragmentLengthModel(options.fragment_mean, options.fragment_st_dev)

        # Set the rng for the fragment length model
        fraglen_model.rng = options.rng

    _LOG.debug("Fragment length model loaded")

    return mut_model, cancer_model, error_model_1, error_model_2, gc_model, fraglen_model


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
        gc_bias_model,
        fraglen_model
    ) = initialize_all_models(options)

    """
    Process Inputs
    """
    _LOG.info(f'Reading {options.reference}.')

    reference_index = SeqIO.index(str(options.reference), "fasta")
    _LOG.debug("Reference file indexed.")

    if _LOG.getEffectiveLevel() < 20:
        count = 0
        for contig in reference_index:
            count += len(reference_index[contig])
        _LOG.debug(f"Length of reference: {count/1_000_000:.2f} Mb")

    input_variants_dict = {x: ContigVariants() for x in reference_index}
    if options.include_vcf:
        _LOG.info(f"Reading input VCF: {options.include_vcf}.")
        if options.cancer:
            sample_names = parse_input_vcf(input_variants_dict,
                                           options.include_vcf,
                                           options.ploidy,
                                           mut_model.homozygous_freq,
                                           reference_index,
                                           options,
                                           tumor_normal=True)

            tumor_ind = sample_names['tumor_sample']
            normal_ind = sample_names['normal_sample']
        else:
            sample_names = parse_input_vcf(input_variants_dict,
                                           options.include_vcf,
                                           options.ploidy,
                                           mut_model.homozygous_freq,
                                           reference_index,
                                           options)

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
    ) = parse_beds(options, reference_index, mut_model.avg_mut_rate)

    if any(bed_files):
        _LOG.debug("Finished reading input beds.")

    # Prepare headers
    bam_header = None
    if options.produce_bam:
        # This is a dictionary that is the list of the contigs and the length of each.
        # This information will be needed later to create the bam header.
        bam_header = {key: len(reference_index[key]) for key in reference_index}

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

    all_variants = {}  # dict of all ContigVariants objects, indexed by contig, which we will collect at the end.
    vcf_files = []
    fasta_files = []
    fastq_files = []

    sam_reads_files = []

    print_fasta_tell = False

    for contig in breaks:

        _LOG.info(f"Generating variants for {contig}")

        # Todo genericize breaks

        input_variants = input_variants_dict[contig]
        # TODO: add the ability to pick up input variants here from previous loop

        local_reference = reference_index[contig]

        # _LOG.info(f'Creating trinucleotide map for {contig}...')
        # local_trinuc_map = map_chromosome(local_reference, mut_model)

        # Since we're only running single threaded for now:
        threadidx = 1

        local_variant_file = options.temp_dir_path / f'{options.output.stem}_tmp_{contig}_{threadidx}.vcf.gz'
        local_fasta_file = None
        if options.produce_fasta:
            local_fasta_file = options.temp_dir_path / f'{options.output.stem}_tmp_{contig}_{threadidx}.fasta'

        _LOG.debug(f'local vcf filename = {local_variant_file}')

        local_bam_pickle_file = None
        if options.produce_bam:
            local_bam_pickle_file = options.temp_dir_path / f'{options.output.stem}_tmp_{contig}_{threadidx}.p.gz'

        if threadidx == 1:
            # init_progress_info()
            pass

        if options.paired_ended:
            max_qual_score = max(seq_error_model_1.quality_scores + seq_error_model_2.quality_scores)
        else:
            max_qual_score = max(seq_error_model_1.quality_scores)

        local_variants = generate_variants(reference=local_reference,
                                           mutation_rate_regions=mutation_rate_dict[contig],
                                           existing_variants=input_variants,
                                           mutation_model=mut_model,
                                           max_qual_score=max_qual_score,
                                           options=options)

        _LOG.info(f'Outputting temp vcf for {contig} for later use')
        if options.produce_fasta:
            _LOG.info(f'Outputting temp fasta for {contig} for later use')
        # This function produces the variant file and the fasta file, if requested
        # TODO pickle dump the ContigVariants object instead. Combine them into one fasta/vcf
        #     at the end.
        write_local_file(local_variant_file,
                         local_variants,
                         local_reference,
                         target_regions_dict[contig],
                         discard_regions_dict[contig],
                         options,
                         local_fasta_file)
        vcf_files.append(local_variant_file)
        fasta_files.append(local_fasta_file)

        if options.produce_fastq or options.produce_bam:
            read1_fastq, read2_fastq = \
                generate_reads(local_reference,
                               local_bam_pickle_file,
                               seq_error_model_1,
                               seq_error_model_2,
                               gc_bias_model,
                               fraglen_model,
                               local_variants,
                               options.temp_dir_path,
                               target_regions_dict[contig],
                               discard_regions_dict[contig],
                               options,
                               contig)

            fastq_files.append((read1_fastq, read2_fastq))
            if options.produce_bam:
                sam_reads_files.append(local_bam_pickle_file)

    if options.produce_vcf:
        _LOG.info(f"Outputting golden vcf: {str(output_file_writer.vcf_fn)}")
        output_file_writer.merge_temp_vcfs(vcf_files)

    if options.produce_fasta:
        _LOG.info(f"Outputting fasta file(s): {', '.join([str(x) for x in output_file_writer.fasta_fns]).strip(', ')}")
        output_file_writer.merge_temp_fastas(fasta_files)

    if options.produce_fastq:
        if options.paired_ended:
            _LOG.info(f"Outputting fastq files: {', '.join([str(x) for x in output_file_writer.fastq_fns]).strip(', ')}")
        else:
            _LOG.info(f"Outputting fastq file: {output_file_writer.fastq_fns[0]}")
        output_file_writer.merge_temp_fastqs(fastq_files, options.paired_ended, options.rng)

    if options.produce_bam:
        _LOG.info(f"Outputting golden bam file: {str(output_file_writer.bam_fn)}")
        contig_list = list(reference_index)
        contigs_by_index = {contig_list[n]: n for n in range(len(contig_list))}
        output_file_writer.output_bam_file(sam_reads_files, contigs_by_index)
