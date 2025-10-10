"""
Runner for generate_reads task
"""
import gzip
import logging
import pickle
import time
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

from pathlib import Path

from Bio import SeqIO

from .utils import Options, parse_input_vcf, parse_beds, OutputFileWriter
from ..common import validate_input_path, validate_output_path
from .single_runner import read_simulator_single
from ..models import MutationModel, SequencingErrorModel, TraditionalQualityModel, FragmentLengthModel
from ..variants import ContigVariants
from .utils.split_inputs import main as split_main
from .utils.stitch_outputs import main as stitch_main

__all__ = ["read_simulator_runner"]

EXTENSIONS = ["gz", "fastq", "bam", "vcf"]

_LOG = logging.getLogger(__name__)

def read_simulator_runner(config: str, output_dir: str, file_prefix: str):
    """
    Run the generate_reads function, which generates simulated mutations in a dataset and corresponding files.

    :param config: This is a configuration file. Keys start with @ symbol. Everything else is ignored.
    :param output_dir: This is the directory where the output will be written.
    :param file_prefix: This is the prefix of the output file names.


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
    analysis_start = time.time()
    _LOG.debug(f'config = {config}')
    _LOG.debug(f'output_dir = {output_dir}')

    _LOG.info(f'Using configuration file {config}')
    config = Path(config).resolve()
    validate_input_path(config)

    # prepare output
    output_dir = Path(output_dir).resolve()
    _LOG.info(f"Saving output files to {output_dir}")

    if not output_dir.is_dir():
        _LOG.info('Creating output dir')
        output_dir.mkdir(parents=True, exist_ok=True)

    # Read options file
    options = Options.from_cli(
        output_dir,
        str(file_prefix),
        config,

    )

    # Set up the recombination later
    reference_index = SeqIO.index(str(options.reference), "fasta")
    reference_keys_with_lens = {key: len(value) for key, value in reference_index.items()}

    count = 0
    for contig in reference_keys_with_lens:
        count += reference_keys_with_lens[contig]
    _LOG.debug(f"Length of total reference: {count / 1_000_000:.2f} Mb")

    # initialize models for run
    (
        mut_model,
        seq_error_model,
        qual_score_model,
        fraglen_model
    ) = initialize_all_models(options)

    input_variants_dict = {x: ContigVariants() for x in reference_keys_with_lens}
    if options.include_vcf:
        _LOG.info(f"Reading input VCF: {options.include_vcf}.")
        parse_input_vcf(
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

    # Prepare headers for final vcf
    bam_header = None
    if options.produce_bam:
        # This is a dictionary that is the list of the contigs and the length of each.
        # This information will be needed later to create the bam header.
        bam_header = reference_keys_with_lens

    if any(bed_files):
        _LOG.debug("Finished reading input beds.")

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf.
    output_file_writer = OutputFileWriter(options=options, bam_header=bam_header)

    # Split file by chunk for parallel analysis or by contig for either parallel or single analysis
    _LOG.info("Splitting reference...")
    splits_files_dict = split_main(options, fraglen_model.fragment_mean, reference_index)

    if options.threads > 1:
        _LOG.info(f"[parallel] Launching {len(splits_files_dict)} NEAT job(s) (max {options.threads} in parallel)...")
    else:
        _LOG.info(f"Launching NEAT!")

    output_opts: list = []
    output_files: list = []
    # a dict with contig keys, then the thread index, and finally the applicable contig variants as the value
    all_variants: dict[str, dict[int, ContigVariants]] = {chrom: {} for chrom in reference_index.keys()}
    for contig in splits_files_dict:
        contig_index = list(reference_keys_with_lens.keys()).index(contig)
        for (thread_idx, splits_file) in enumerate(splits_files_dict[contig]):
            current_output_dir = options.temp_dir_path / splits_file.stem
            current_output_dir.mkdir(parents=True, exist_ok=True)
            # Create local filenames based on fasta indexing scheme.
            fq1 = None
            fq2 = None
            bam = None
            vcf = None
            if options.produce_fastq:
                fq1 = current_output_dir / options.fq1.name
                # validate to double-check we don't have name collisions
                validate_output_path(fq1, True, False)
                if options.paired_ended:
                    # if paired ended, there must be a fq2
                    fq2 = current_output_dir / options.fq2.name
                    validate_output_path(fq2, True, False)
            if options.produce_bam:
                bam = current_output_dir / options.bam.name
                validate_output_path(bam, True, False)
            if options.produce_vcf:
                vcf = current_output_dir / options.vcf.name
                validate_output_path(vcf, True, False)

            current_options = options.copy_with_changes(splits_file, current_output_dir, fq1, fq2, vcf, bam)
            if options.threads == 1:
                idx, contig, files_written, local_variants = read_simulator_single(
                    1,
                    current_options,
                    contig,
                    contig_index,
                    mut_model,
                    seq_error_model,
                    qual_score_model,
                    fraglen_model,
                    input_variants_dict[contig],
                    target_regions_dict[contig],
                    discard_regions_dict[contig],
                    mutation_rate_dict[contig],
                )
                _LOG.info(f"Completed simulating contig {contig}.")
                all_variants[contig][idx] = local_variants
                output_files.append((thread_idx, files_written))
            else:
                output_opts.append((
                    thread_idx,
                    current_options,
                    contig,
                    contig_index,
                    mut_model,
                    seq_error_model,
                    qual_score_model,
                    fraglen_model,
                    input_variants_dict[contig],
                    target_regions_dict[contig],
                    discard_regions_dict[contig],
                    mutation_rate_dict[contig],
                ))

    if options.threads > 1:
        pool = ThreadPool(options.threads)
        (thread_idx, contig_name, files_written, local_variants) = pool.starmap(read_simulator_single, output_opts)
        all_variants[contig_name][thread_idx] = local_variants
        output_files.append((thread_idx, files_written))

    # stitch fastq and bams
    stitch_main(options, output_files)

    # sort all variants and write out final VCF
    if options.produce_vcf:
        write_final_vcf(all_variants, reference_index, output_file_writer)

    _LOG.info(f"Read simulator complete in {time.time() - analysis_start} s")

def write_final_vcf(
        all_variants: dict[str, dict[int, ContigVariants]],
        ref_index: dict,
        ofw: OutputFileWriter,
):
    # take contig name from ref index because we know it is in the proper order
    for contig in ref_index.keys():
        block_vars_dict = all_variants[contig]
        # Ensuring these are put together in the right order, given that they might have
        #    been written out of order, due to race conditions.
        sorted_keys = sorted(block_vars_dict.keys())
        for block in sorted_keys:
            local_variants = block_vars_dict[block]
            locations = sorted(local_variants.variant_locations)
            for location in locations:
                for variant in local_variants[location]:
                    ref, alt = local_variants.get_ref_alt(variant, ref_index[contig])
                    sample = local_variants.get_sample_info(variant)
                    # +1 to position because the VCF uses 1-based coordinates
                    #          .id should give the more complete name
                    line = f"{ref_index[contig].id}\t" \
                           f"{variant.position1 + 1}\t" \
                           f"{local_variants.generate_field(variant, 'ID')}\t" \
                           f"{ref}\t" \
                           f"{alt}\t" \
                           f"{variant.qual_score}\t" \
                           f"{local_variants.generate_field(variant, 'FILTER')}\t" \
                           f"{local_variants.generate_field(variant, 'INFO')}\t" \
                           f"{local_variants.generate_field(variant, 'FORMAT')}\t" \
                           f"{sample}\n"
                    ofw.write_vcf_record(line)

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

    _LOG.debug("Mutation models loaded")

    # We need sequencing errors to get the quality score attributes, even for the vcf
    if options.error_model:
        error_models = pickle.load(gzip.open(options.error_model))
        error_model = error_models["error_model1"]
        quality_score_model = error_models["qual_score_model1"]
    else:
        # Use all the default values
        error_model = SequencingErrorModel()
        quality_score_model = TraditionalQualityModel()

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
        error_model, \
        quality_score_model, \
        fraglen_model