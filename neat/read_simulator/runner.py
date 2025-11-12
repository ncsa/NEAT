"""
Runner for generate_reads task
"""
import gzip
import logging
import os
import pickle
import subprocess
import time
import multiprocessing as mp
from math import ceil

from pathlib import Path

import pysam
from pysam import bcftools
from Bio import SeqIO

from .utils import Options, OutputFileWriter, parse_beds, parse_input_vcf
from ..common import validate_input_path, validate_output_path
from .single_runner import read_simulator_single
from ..models import SequencingErrorModel
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

    # We need sequencing errors to get the quality score attributes, even for the vcf
    if options.avg_seq_error:
        average_error = options.avg_seq_error
    elif options.error_model:
        error_models = pickle.load(gzip.open(options.error_model))
        average_error = error_models["error_model1"].average_error
        # We just need the error value
        del error_models
    else:
        # Use the default value
        average_error = 0.009228843915252066

    # _LOG.debug('Sequencing error and quality score models loaded')
    # We need to estimate how many total errors to add
    total_reference_length = sum(reference_keys_with_lens.values())
    # This is an estimate due to some random effects
    number_of_bases_in_analysis = total_reference_length * options.coverage
    # Each base called has an "average_error" chance of being an error
    total_errors = ceil(average_error * number_of_bases_in_analysis)
    # Normalization gives the percent value for each item with a sum of all values being 1.0
    normalized_counts = {k: v / total_reference_length for (k, v) in reference_keys_with_lens.items()}
    # Multiply the normalized count by the total errors. This gives the number of errors that should be
    # introduced into the contig
    errors_per_contig = {k: ceil(v * total_errors) for (k, v) in normalized_counts.items()}

    count = 0
    for contig in reference_keys_with_lens:
        count += reference_keys_with_lens[contig]
    _LOG.debug(f"Length of total reference: {count / 1_000_000:.2f} Mb")

    # Note that parse_beds will return None for any empty or non-input files
    # parse input targeted regions, if present.
    if any((options.target_bed, options.discard_bed, options.mutation_bed)):
        _LOG.info(f"Reading input bed files.")

    (
        target_regions_dict,
        discard_regions_dict,
        mutation_rate_dict
    ) = parse_beds(options, reference_keys_with_lens)

    input_variants_dict = {x: ContigVariants() for x in reference_keys_with_lens}
    if options.include_vcf:
        _LOG.info(f"Reading input VCF: {options.include_vcf}.")
        parse_input_vcf(
            input_variants_dict,
            options.include_vcf,
            options.ploidy,
            reference_index,
            options
        )

        _LOG.debug("Finished reading input vcf file")

    if any((options.target_bed, options.discard_bed, options.mutation_bed)):
        _LOG.debug("Finished reading input beds.")

    # Prepare headers for final vcf
    bam_header = None
    if options.produce_bam:
        # This is a dictionary that is the list of the contigs and the length of each.
        # This information will be needed later to create the bam header.
        bam_header = reference_keys_with_lens

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf. We create the overall bam with no header, as it will get a header from
    # merging the smaller bams.
    output_file_writer = OutputFileWriter(options=options, vcf_header = reference_keys_with_lens, bam_header=None)

    # Split file by chunk for parallel analysis or by contig for either parallel or single analysis
    _LOG.info("Splitting reference...")
    (splits_files_dict, count) = split_main(
        options,
        reference_index
    )

    if options.threads > 1:
        _LOG.info(f"[parallel] Launching {count} NEAT job(s) (max {options.threads} in parallel)...")
    else:
        _LOG.info(f"Launching NEAT!")

    output_opts: list = []
    output_files: list = []
    # a dict with contig keys, then the thread index, and finally the applicable contig variants as the value
    # TODO Remove if not needed
    # all_variants: dict[str, dict[int, ContigVariants]] = {chrom: {} for chrom in reference_index.keys()}
    thread_idx = 1
    contig_list = list(reference_keys_with_lens.keys())
    contig_dict = {contig: contig_list.index(contig) for contig in reference_keys_with_lens.keys()}
    for contig in splits_files_dict:
        contig_index = contig_dict[contig]
        for ((start, length), splits_file) in splits_files_dict[contig].items():
            block_percentage = length / reference_keys_with_lens[contig]
            block_errors = errors_per_contig[contig] * block_percentage
            estimated_number_of_reads = (length // options.read_len) * options.coverage
            errors_per_read = round(block_errors / estimated_number_of_reads)
            if errors_per_read < 1.0 and block_errors > 0:
                # We know we need a few errors, but it's a small number total
                if options.rng.random() < average_error:
                    errors_per_read += 1
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
            current_options = options.copy_with_changes(splits_file, current_output_dir, fq1, fq2, bam, vcf)
            if options.threads == 1:
                idx, contig, local_variants, files_written = read_simulator_single(
                    1,
                    start,
                    current_options,
                    bam_header,
                    contig,
                    contig_index,
                    input_variants_dict[contig],
                    target_regions_dict[contig],
                    discard_regions_dict[contig],
                    mutation_rate_dict[contig],
                    errors_per_read,
                )
                _LOG.info(f"Completed simulating contig {contig}.")
                # TODO Remove if not needed
                # all_variants[contig][idx] = local_variants
                output_files.append((thread_idx, files_written))
            else:
                thread_input_variants = filter_thread_variants(input_variants_dict[contig], (start, start+length))
                thread_target_regions = filter_bed_regions(target_regions_dict[contig], (start, start+length))
                thread_discard_regions = filter_bed_regions(discard_regions_dict[contig], (start, start + length))
                thread_mutation_regions = filter_bed_regions(mutation_rate_dict[contig], (start, start + length))

                output_opts.append((
                    thread_idx,
                    start,
                    current_options,
                    bam_header,
                    contig,
                    contig_index,
                    thread_input_variants,
                    thread_target_regions,
                    thread_discard_regions,
                    thread_mutation_regions,
                    errors_per_read,
                ))
            thread_idx += 1

    if options.threads > 1:
        pool = mp.Pool(options.threads)
        results = pool.starmap_async(read_simulator_single, output_opts)

        _LOG.info(f"launching mutiprocess simulation, recording results.")
        pool.close()
        pool.join()

        # Need to organize the results, as above
        for idx, contig, local_variants, files_written in results.get():
            # TODO Remove if not needed
            # all_variants[contig][idx] = local_variants
            output_files.append((idx, files_written))

    _LOG.info("Processing complete, writing output")

    start = time.time()
    if options.produce_bam:
        stitch_main(output_file_writer, output_files, options.threads)
        _LOG.info(f"It took {time.time() - start} s to write the bam file")
    else:
        stitch_main(output_file_writer, output_files)

    output_file_writer.flush_and_close_files()
    force = False
    if options.overwrite_output:
        force = True
    for file in options.output_files:
        if file.suffix == ".bam":
            _LOG.info(f"bam file: {file}")
            pysam.index(f"{str(file)}", "-@", str(options.threads))
            _LOG.info(f'bam index: {file}.bai')
        elif "fastq" in file.name:
            continue
        else:
            _LOG.info("Sorting VCF file")
            temp_file = str(options.temp_dir_path / "temp.sorted.vcf.gz")
            subprocess.run(["bcftools", "sort", "-o", temp_file, "-Ob9", str(file)])
            Path(temp_file).is_file()
            os.rename(temp_file, str(file))
            _LOG.info("Indexing vcf")
            pysam.tabix_index(str(file), preset="vcf", force=force)

    _LOG.info(f"Read simulator complete in {time.time() - analysis_start} s")

def filter_thread_variants(contig_variants: ContigVariants, coords: tuple[int, int]) -> ContigVariants:
    ret_contig_vars = ContigVariants()
    for variant_loc in contig_variants.variant_locations:
        if coords[0] <= variant_loc < coords[1]:
            # The variant location is within our area of interest
            variants_of_interest = contig_variants[variant_loc]
            for variant in variants_of_interest:
                ret_contig_vars.add_variant(variant)

    return ret_contig_vars

def filter_bed_regions(regions: list, coords: tuple[int,int]) -> list:
    ret_list = []
    for region in regions:
        if region[0] <= coords[0] < region[1] or \
            region[0] < coords[1] <= region[1] or \
            coords[0] <= region[0] < coords[1] or \
            coords[0] < region[1] <= coords[1]:
            ret_list.append(region)
    return ret_list