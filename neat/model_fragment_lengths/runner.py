"""
Creates a model of fragment lengths in a dataset
"""

import gzip
import pickle
import numpy as np
import logging
import sys

from pathlib import Path

from .utils import count_frags, filter_lengths
from ..models import FragmentLengthModel
from ..common import validate_output_path

__all__ = [
    "compute_fraglen_runner"
]

_LOG = logging.getLogger(__name__)


def compute_fraglen_runner(
        file: str | Path,
        filter_minreads: int,
        output_dir: str | Path,
        output_prefix: str,
        overwrite: bool = False):
    """
    Main function takes 2 arguments:
    :param file: a path to a sam or bam file_list input
    :param filter_minreads: minimum number of reads needed to count a fragment length.
                            If 0, then filtering will be skipped
    :param output_dir: The directory where to write files
    :param output_prefix: the prefix for filenames
    :param overwrite: If true, any existing file_list found will be overwritten
    """
    _LOG.info("Generating fragment length model")

    input_file = file
    out_file_name = output_prefix + ".pickle.gz"
    output_file = Path(output_dir) / out_file_name
    validate_output_path(output_file, True, overwrite)

    _LOG.debug(f"File: {input_file}")
    _LOG.debug(f'Minimum reads: {filter_minreads}')
    _LOG.debug(f'Output: {output_file}')
    _LOG.debug(f'Overwrite: {overwrite}')

    _LOG.info("Counting fragments")
    all_tlens = count_frags(input_file)
    if not all_tlens:
        _LOG.error("No valid template lengths in sam file_list.")
        sys.exit(1)

    _LOG.info("Filtering fragments")
    filtered_lengths = filter_lengths(all_tlens, filter_minreads)

    if not filtered_lengths:
        _LOG.error("No data passed the filter, nothing to calculate. Try adjusting the filter settings.")
        sys.exit(1)

    _LOG.info("Building model")
    st_dev = float(np.std(filtered_lengths))
    mean = float(np.mean(filtered_lengths))

    model = FragmentLengthModel(mean, st_dev)
    _LOG.info(f'Saving model: {output_file}')
    with gzip.open(output_file, 'w+') as outfile:
        pickle.dump(model, outfile)

    _LOG.info("Modeling complete.")
