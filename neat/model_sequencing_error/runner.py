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
        files: list,
        offset: int,
        qual_scores: list | int,
        max_reads: int,
        overwrite: bool,
        output_prefix: str,
        pileup: str = None,
        plot: bool = False,
):
    """
    Sets up and calls the sequencing error modeling core code.

    :param files: This is the input file or pair of input files, which can be in fastq (bgzipped or not) format.
    :param offset: This is the quality offset. Illumina machines use the sanger model qual score + 33 to
        derive the character for the quality array.
    :param qual_scores: These are the possible quality scores for the dataset. This can either be a list or a single
        integer.
    :param max_reads: The maximum number of reads to process. This can speed up the time taken to create the model,
        at the expense of accuracy.
    :param overwrite: True to overwrite input, mainly a debugging option
    :param output_prefix: The name of the file to write the output.
    :param pileup: If pileup file included, get stats from that
    :param plot: run optional plotting.
    """

    if len(files) > 2:
        _LOG.info(f'Only processing first two input files')
        files = files[:2]

    for file in files:
        validate_input_path(file)

    _LOG.debug(f'Input Files: {", ".join([str(x) for x in files])}')

    _LOG.debug(f"Quality offset: {offset}")

    final_quality_scores: list
    if len(qual_scores) == 1:
        final_quality_scores = list(range(1, qual_scores[0] + 1))
    else:
        final_quality_scores = sorted(qual_scores)

    _LOG.debug(f'Quality scores: {final_quality_scores}')
    if max_reads == -1:
        num_records_to_process = np.inf
    else:
        num_records_to_process = max_reads

    _LOG.debug(f'Maximum number of records to process: '
               f'{"all" if num_records_to_process == np.inf else num_records_to_process}')

    if pileup:
        pileup = Path(pileup)
        validate_input_path(pileup)

    _LOG.debug(f'Pileup file: {pileup}')
    _LOG.debug(f"Plot the data? {plot}")
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
    for file in files:
        file_num += 1
        _LOG.info(f'Reading file {file_num} of {len(files)}')
        parameters_by_position, file_avg_error, read_length = parse_file(file,
                                                                         final_quality_scores,
                                                                         num_records_to_process,
                                                                         offset,
                                                                         read_length)

        read_parameters.append(parameters_by_position)
        average_errors.append(file_avg_error)

        _LOG.info(f'Finished reading file {file_num}')

    read_parameters = np.asarray(read_parameters)
    average_error = np.average(average_errors)

    _LOG.info(f"Found an average error of {average_error} across {len(files)} file(s).")

    if plot:
        _LOG.info("Plotting coming soon! Sorry!")

    if pileup:
        _LOG.info("Pileup features not yet available. Using default parameters")

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

    # First model, always produced
    seq_err_model = SequencingErrorModel(
        avg_seq_error=average_error,
        read_length=read_length,
        quality_scores=np.array(final_quality_scores),
        qual_score_probs=read_parameters[0],
    )

    if len(files) == 1:
        with gzip.open(output_file, 'w') as out_model:
            pickle.dump([seq_err_model, None], out_model)

    else:
        # Second model if a second input was given
        seq_err_model_r2 = SequencingErrorModel(
            avg_seq_error=average_error,
            read_length=read_length,
            quality_scores=np.array(final_quality_scores),
            qual_score_probs=read_parameters[1]
        )
        with gzip.open(output_file, 'w') as out_model:
            pickle.dump([seq_err_model, seq_err_model_r2], out_model)

    _LOG.info("Modeling sequencing errors is complete, have a nice day.")
