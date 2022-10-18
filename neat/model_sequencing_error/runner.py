"""
Creates a model of sequencing errors in the dataset. One assumption of this model is that low quality scores indicate
sequencing errors. We should verify this is true in practice and/or how that assumption lines up with variant callers.

I wonder if it would be possible to get a list of all items a variant caller deems errors.
"""

import gzip
import pickle
import numpy as np
import logging

from pathlib import Path

from .utils import parse_file
from ..common import validate_output_path, validate_input_path
from ..models import SequencingErrorModel

__all__ = [
    "model_seq_err_runner"
]

_LOG = logging.getLogger(__name__)


def model_seq_err_runner(
        file_list: list,
        offset: int,
        qual_scores: list | int,
        max_reads: int,
        pileup: str,
        plot: bool,
        overwrite: bool,
        output_prefix: str
):
    """
    Sets up and calls the sequencing error modeling core code.

    :param file_list: This is the input file, which can be in fastq (gzipped or not), sam, or bam format.
    :param offset: This is the quality offset. Illumina machines use the sanger model qual score + 33 to
        derive the character for the quality array.
    :param qual_scores: These are the possible quality scores for the dataset. This can either be a list or a single
        integer.
    :param max_reads: The maximum number or reads to process. This can speed up the time taken to create the model,
        at the expense of accuracy.
    :param pileup: If pileup file included, get stats from that
    :param plot: run optional plotting.
    :param overwrite: True to overwrite input, mainly a debugging option
    :param output_prefix: The name of the file to write the output.
    """

    input_files = []
    for file in file_list:
        input_files.append(Path(file))
    common_fastq_suffixes = ['.fq', '.fastq']
    filetype = None
    suffixes = []
    for file in input_files:
        suffixes.extend(file.suffixes)
    fastq_check = any([x for x in suffixes if x in common_fastq_suffixes])
    bam_check = '.bam' in suffixes
    sam_check = '.sam' in suffixes
    if sam_check or bam_check:
        # checking for mixing of file types.
        if fastq_check:
            raise ValueError("NEAT detected mixed file types for the input. Please use consistent file types.")
        elif sam_check or bam_check:
            filetype = 'bam'
    elif fastq_check:
        filetype = "fastq"
    else:
        raise ValueError("Unrecognized file extension in input files. Please use .fq or .fastq for fastq files and"
                         ".sam or .bam for sam/bam formatted files.")

    _LOG.debug(f'Input Files: {", ".join([str(x) for x in input_files])}')
    _LOG.debug(f'Filetype: {filetype}')

    _LOG.debug(f"Quality offset: {offset}")

    final_quality_scores: list
    if len(qual_scores) == 1:
        final_quality_scores = list(range(qual_scores[0] + 1))
    else:
        final_quality_scores = qual_scores

    _LOG.debug(f'Quality scores: {final_quality_scores}')
    if max_reads == -1:
        num_records_to_process = np.inf
    else:
        num_records_to_process = max_reads

    _LOG.debug(f'Maximum number of records to process: {num_records_to_process}')

    if pileup:
        pileup = Path(pileup)
        validate_input_path(pileup)

    _LOG.debug(f'Pileup file: {pileup}')
    _LOG.debug(f"Plot the data? {plot}")
    _LOG.debug(f'Overwrite existing data? {overwrite}')

    validate_output_path(output_prefix, is_file=False)
    output_file = Path(output_prefix).with_suffix(".p.gz")
    validate_output_path(output_file, overwrite=overwrite)
    _LOG.info(f'Writing output to: {output_file}')

    probabilities = []
    average_errors = []
    overall_read_length = 0
    for file in input_files:
        file_qual_score_probs, file_avg_error, file_readlen = parse_file(file, filetype, final_quality_scores, offset, num_records_to_process)
        probabilities.append(file_qual_score_probs)
        average_errors.append(file_avg_error)
        if not overall_read_length:
            overall_read_length = file_readlen
        elif file_readlen != overall_read_length:
            _LOG.warning("The input files seem to have different read lengths. Using the shorter of the two")
            overall_read_length = min(overall_read_length, file_readlen)

    final_probs = np.zeros((overall_read_length, len(qual_scores)), dtype=float)
    for i in range(overall_read_length):
        probs_at_location = [x[i] for x in probabilities]
        temp_probs = np.average(probs_at_location, axis=0)
        # technically, temp_probs should sum to 1, but just in case
        final_probs[i] = temp_probs/sum(temp_probs)

    average_error = sum(average_errors)/len(average_errors)
    _LOG.info(f"Found an average error of {average_error} across {len(input_files)} file(s).")

    if plot:
        _LOG.info("Plotting coming soon! Sorry!")

    if pileup:
        _LOG.info("Pileup features not yet available. Using default parameters")

    # Generate the model
    seq_err_model = SequencingErrorModel(
        avg_seq_error=average_error,
        read_length=overall_read_length,
        quality_scores=np.array(qual_scores),
        qual_score_probs=final_probs,
    )
    # finally, let's save our output model
    _LOG.info(f'Saving model: {output_file}')

    pickle.dump(seq_err_model, gzip.open(output_file, 'w'))

    _LOG.info("Modeling sequencing errors is complete, have a nice day.")
