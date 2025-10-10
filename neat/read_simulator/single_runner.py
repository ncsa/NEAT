"""
Runner for read-simulator in single-ended mode
"""
from turtledemo.chaos import coosys

from Bio import SeqIO
import logging
import pickle
import gzip

from .utils import parse_input_vcf, parse_beds, OutputFileWriter, \
    generate_variants, generate_reads, Options
from ..variants import ContigVariants

from ..models import MutationModel, SequencingErrorModel, FragmentLengthModel, TraditionalQualityModel

__all__ = ["read_simulator_single"]

_LOG = logging.getLogger(__name__)

def read_simulator_single(
        thread_idx: int,
        local_options: Options,
        contig_name: str,
        mut_model: MutationModel,
        seq_error_model: SequencingErrorModel,
        qual_score_model: TraditionalQualityModel,
        fraglen_model: FragmentLengthModel,
        input_variants: ContigVariants,
        target_regions: list,
        discard_regions: list,
        mutation_rate_regions: list,
) -> tuple[]:
    """
    inputs:
    :param thread_idx: index of current thread
    :param local_options: options for current thread and reference chunk
    :param contig_name: The original list of contig names.
    :param mut_model: The mutation model for the run
    :param seq_error_model: Sequencing error model for run
    :param qual_score_model: Quality score model for run
    :param fraglen_model: Fragment length model for the run
    :param input_variants: Input variants for the current block
    :param target_regions: Regions to target in analysis
    :param discard_regions: Regions to discard from analysis
    :param mutation_rate_regions: Mutation rate bed regions for this section

    Ideally this should work for either a file chunk or contig. We'll assume here that we're
    getting either an entire contig or a file chunk, and that no new subdivisions are needed.
    We can just read in the file.

    Read input models or default models, as specified by user.
    """
    if thread_idx != 1:
        _THREAD_LOG = logging.getLogger(f"thread_{thread_idx}")
        _THREAD_LOG.propagate = False
    else:
        _THREAD_LOG = _LOG
    """
    Process Inputs
    """
    _THREAD_LOG.info(f'Reading {local_options.reference}.')
    local_ref_index = SeqIO.index(str(local_options.reference), "fasta")

    # For the local bam, we will forgo the header, and then add it at the end.
    bam_header = None

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf.
    # We'll also keep track here of what files we are producing.
    # We don't really need to write out the VCF. We should be able to store it in memory
    local_options.produce_vcf = False
    local_output_file_writer = OutputFileWriter(options=local_options, bam_header=bam_header)
    """
    Begin Analysis
    """
    max_qual_score = max(qual_score_model.quality_scores)

    local_variants = generate_variants(
        reference=local_ref_index,
        mutation_rate_regions=mutation_rate_regions,
        existing_variants=input_variants,
        mutation_model=mut_model,
        max_qual_score=max_qual_score,
        options=local_options,
    )

    if local_options.produce_fastq or local_options.produce_bam:
        generate_reads(
            thread_idx,
            local_ref_index[contig_name],
            seq_error_model,
            qual_score_model,
            fraglen_model,
            local_variants,
            target_regions,
            discard_regions,
            local_options,
            contig_name,
            local_output_file_writer,
        )

    local_output_file_writer.close_files()
    return (
        thread_idx,
        contig_name,
        local_output_file_writer.files_to_write.keys(),
        local_variants
    )
