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
from ..variants import Insertion, Deletion, SingleNucleotideVariant

__all__ = [
    "model_seq_err_runner"
]

_LOG = logging.getLogger(__name__)


def model_seq_err_runner(
        file1: str,
        file2: str,
        offset: int,
        qual_scores: list | int,
        max_reads: int,
        # pileup: str,
        # plot: bool,
        overwrite: bool,
        output_prefix: str
):
    """
    Sets up and calls the sequencing error modeling core code.

    :param file1: This is the input file, which can be in fastq (bgzipped or not) format.
    :param file2: This is a second input file, for paired ended reads, which can be in fastq (bgzipped or not) format.
    :param offset: This is the quality offset. Illumina machines use the sanger model qual score + 33 to
        derive the character for the quality array.
    :param qual_scores: These are the possible quality scores for the dataset. This can either be a list or a single
        integer.
    :param max_reads: The maximum number of reads to process. This can speed up the time taken to create the model,
        at the expense of accuracy.
    :param pileup: If pileup file included, get stats from that
    :param plot: run optional plotting.
    :param overwrite: True to overwrite input, mainly a debugging option
    :param output_prefix: The name of the file to write the output.
    """

    validate_input_path(file1)
    input_files = [file1]
    if file2:
        validate_input_path(file2)
        input_files.append(file2)

    _LOG.debug(f'Input Files: {", ".join([str(x) for x in input_files])}')

    _LOG.debug(f"Quality offset: {offset}")

    final_quality_scores: list
    if len(qual_scores) == 1:
        final_quality_scores = list(range(1, qual_scores[0] + 1))
    else:
        final_quality_scores = qual_scores

    _LOG.debug(f'Quality scores: {final_quality_scores}')
    if max_reads == -1:
        num_records_to_process = np.inf
    else:
        num_records_to_process = max_reads

    _LOG.debug(f'Maximum number of records to process: {num_records_to_process}')

    # if pileup:
    #     pileup = Path(pileup)
    #     validate_input_path(pileup)
    #
    # _LOG.debug(f'Pileup file: {pileup}')
    # _LOG.debug(f"Plot the data? {plot}")
    _LOG.debug(f'Overwrite existing data? {overwrite}')

    validate_output_path(output_prefix, is_file=False)
    output_prefix = Path(output_prefix)

    """
    Previously, this was
    
    output_file = Path(output_prefix).with_suffix('.p.gz')
    
    But this tended to drop parts of the name, if they included dots: 
    e.g., output_prefix = "/path/to/samplename.testrun" outputted to /path/to/samplename.p.gz
    
    This way avoids that issue and outputs to /path/to/samplename.testrun.p.gz, 
    even though it is a little more cumbersome
    """
    output_file = output_prefix.parent / f'{output_prefix.name}.p.gz'
    validate_output_path(output_file, overwrite=overwrite)
    _LOG.info(f'Writing output to: {output_file}')

    read_parameters = []
    average_errors = []
    read_length = 0
    file_num = 0
    for file in input_files:
        file_num += 1
        _LOG.info(f'Reading file {file_num} of {len(input_files)}')
        parameters_by_position, file_avg_error, file_readlen = parse_file(file,
                                                                          final_quality_scores,
                                                                          num_records_to_process,
                                                                          offset)
        read_parameters.append(parameters_by_position)
        average_errors.append(file_avg_error)
        if not read_length:
            read_length = file_readlen
        elif file_readlen != read_length:
            _LOG.warning("Read lengths inconsistent between reads. Using the smaller value for the model")
            read_length = min(file_readlen, read_length)

        _LOG.info(f'Finished reading file {file_num}')

    read_parameters = np.asarray(read_parameters)
    average_error = np.average(average_errors)

    _LOG.info(f"Found an average error of {average_error} across {len(input_files)} file(s).")

    # if plot:
    #     _LOG.info("Plotting coming soon! Sorry!")
    #
    # if pileup:
    #     _LOG.info("Pileup features not yet available. Using default parameters")

    # Generate and save the model

    # Default values from the original NEAT
    # TODO incorporate these into the calculations
    error_transition_matrix = np.array(
        [[0.0, 0.4918, 0.3377, 0.1705],
         [0.5238, 0.0, 0.2661, 0.2101],
         [0.3755, 0.2355, 0.0, 0.389],
         [0.2505, 0.2552, 0.4943, 0.0]]
    )

    error_variant_probs = {Insertion: 0.004, Deletion: 0.006, SingleNucleotideVariant: 0.99}
    indel_len_model = {1: 0.999, 2: 0.001}
    insertion_model = np.array([0.25, 0.25, 0.25, 0.25])

    _LOG.info(f'Saving model: {output_file}')
    with gzip.open(output_file, 'w') as out_model:
        pickle.dump(
            SequencingErrorModel(
                avg_seq_error=average_error,
                read_length=read_length,
                transition_matrix=error_transition_matrix,
                quality_scores=np.array(final_quality_scores),
                qual_score_probs=read_parameters[0],
                variant_probs=error_variant_probs,
                indel_len_model=indel_len_model,
                insertion_model=insertion_model
            ),
            out_model
        )

        if len(input_files) == 2:
            pickle.dump(
                SequencingErrorModel(
                    avg_seq_error=average_error,
                    read_length=read_length,
                    transition_matrix=error_transition_matrix,
                    quality_scores=np.array(final_quality_scores),
                    qual_score_probs=read_parameters[1],
                    variant_probs=error_variant_probs,
                    indel_len_model=indel_len_model,
                    insertion_model=insertion_model
                ),
                out_model
            )

    _LOG.info("Modeling sequencing errors is complete, have a nice day.")
