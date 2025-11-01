"""
Runner for read-simulator in single-ended mode
"""
import gzip
import os
import pickle

import pysam
from Bio import SeqIO, bgzf
import logging
from pathlib import Path

from .utils import OutputFileWriter, \
    generate_variants, generate_reads, Options, recalibrate_mutation_regions
from ..variants import ContigVariants

from ..models import MutationModel, SequencingErrorModel, FragmentLengthModel, TraditionalQualityModel

__all__ = ["read_simulator_single"]

_LOG = logging.getLogger(__name__)

def read_simulator_single(
        thread_idx: int,
        block_start: int,
        local_options: Options,
        bam_header: list | None,
        contig_name: str,
        contig_index: int,
        input_variants_local: ContigVariants,
        target_regions: list,
        discard_regions: list,
        mutation_regions: list,
) -> tuple[int, str, ContigVariants, dict[str, Path], ]:
    """
    inputs:
    :param thread_idx: index of current thread
    :param block_start: Where on the reference does this block start? For a full contig, this will be 0.
    :param local_options: options for current thread and reference chunk
    :param bam_header: Pass in the outer bam header.
    :param contig_name: The original list of contig names.
    :param contig_index: The index of the contig which this chunk comes from
    :param input_variants_local: The input variants for this block
    TODO I'm counting on the target and discard regions not being used. They are likely broken with multithreading.
            Probably they make more sense after the fact now
    :param target_regions: Target regions for the run
    :param discard_regions:  discard regions for the run
    :param mutation_regions: mutation regions (unchecked) for the run

    Ideally this should work for either a file chunk or contig. We'll assume here that we're
    getting either an entire contig or a file chunk, and that no new subdivisions are needed.
    We can just read in the file.

    Read input models or default models, as specified by user.
    """
    # if thread_idx != 1:
    #     _THREAD_LOG = logging.getLogger(f"thread_{thread_idx}")
    #     _THREAD_LOG.propagate = False
    # else:
    #     _THREAD_LOG = _LOG

    """
    Load models (note that this means each thread has it's very own copy of the models. It also means they may be trying 
    to read the initial beds at the same time. But not sure of a better solution at this point.
    """
    # _THREAD_LOG.info('Initializing models')
    # initialize models for run
    (
        mut_model,
        seq_error_model,
        qual_score_model,
        fraglen_model
    ) = initialize_all_models(local_options)

    """
    Process Inputs
    """
    # _THREAD_LOG.info(f'Reading {local_options.reference}.')
    local_ref_index = SeqIO.index(str(local_options.reference), "fasta")
    local_ref_name = list(local_ref_index.keys())[0]
    local_seq_record = local_ref_index[local_ref_name]

    if len(local_seq_record) < local_options.read_length:
        _LOG.debug("Record too small for processing")


    coords = (block_start, block_start+len(local_seq_record))
    mutation_rate_regions = recalibrate_mutation_regions(mutation_regions, coords, mut_model.avg_mut_rate)

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf.
    # We'll also keep track here of what files we are producing.
    # We don't really need to write out the VCF. We should be able to store it in memory
    local_options.produce_vcf = False
    local_output_file_writer = OutputFileWriter(options=local_options, header=bam_header)
    """
    Begin Analysis
    """
    max_qual_score = max(qual_score_model.quality_scores)

    local_variants = generate_variants(
        reference=local_seq_record,
        ref_start=block_start,
        mutation_rate_regions=mutation_rate_regions,
        existing_variants=input_variants_local,
        mutation_model=mut_model,
        max_qual_score=max_qual_score,
        options=local_options,
    )

    if local_options.produce_fastq or local_options.produce_bam:
        reads_to_write = generate_reads(
            thread_idx,
            local_seq_record,
            seq_error_model,
            qual_score_model,
            fraglen_model,
            local_variants,
            target_regions,
            discard_regions,
            local_options,
            local_output_file_writer,
            contig_name,
            contig_index,
        )
        if local_options.produce_bam:
            # Writing an intermediate bam that is sorted, to make compiling them together at the end easier.
            bam_handle = local_output_file_writer.files_to_write[local_output_file_writer.bam]
            for read_data in reads_to_write:
                read1 = read_data[0]
                read2 = read_data[1]
                if read1:
                    local_output_file_writer.write_bam_record(
                        read1,
                        contig_index,
                        bam_handle,
                        local_options.read_len
                    )
                if read2:
                    local_output_file_writer.write_bam_record(
                        read2,
                        contig_index,
                        bam_handle,
                        local_options.read_len
                    )
            bam_handle.flush()
            bam_handle.close()
            sorted_bam = local_output_file_writer.bam.with_suffix(".sorted.bam")
            pysam.sort("-@", str(local_options.threads), "-o", str(sorted_bam), str(local_output_file_writer.bam))
            os.rename(str(sorted_bam), str(local_output_file_writer.bam))
            _LOG.info(f"bam for thread {thread_idx} written")

    write_block_vcf(local_variants, contig_name, local_ref_index, local_output_file_writer)

    local_output_file_writer.flush_and_close_files(False)
    file_dict = {
        "fq1": local_output_file_writer.fq1,
        "fq2": local_output_file_writer.fq2,
        "vcf": local_output_file_writer.vcf,
        "bam": local_options.bam,
    }
    return (
        thread_idx,
        contig_name,
        local_variants,
        file_dict,
    )

def write_block_vcf(
        local_variants: ContigVariants,
        contig: str,
        ref_index: dict,
        ofw: OutputFileWriter,
):
    # take contig name from ref index because we know it is in the proper order
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

    # _LOG.debug("Mutation models loaded")

    # We need sequencing errors to get the quality score attributes, even for the vcf
    if options.error_model:
        error_models = pickle.load(gzip.open(options.error_model))
        error_model = error_models["error_model1"]
        quality_score_model = error_models["qual_score_model1"]
    else:
        # Use all the default values
        error_model = SequencingErrorModel()
        quality_score_model = TraditionalQualityModel()

    # _LOG.debug('Sequencing error and quality score models loaded')

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

    # _LOG.debug("Fragment length model loaded")

    return \
        mut_model, \
        error_model, \
        quality_score_model, \
        fraglen_model
