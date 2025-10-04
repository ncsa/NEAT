"""
Runner for generate_reads task
"""
import logging

from pathlib import Path

from .utils import Options
from ..common import validate_input_path
from .parallel_runner import main as parallel_runner
from .single_runner import read_simulator_single as single_runner

__all__ = ["read_simulator_runner"]

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

    # set the output dir plus file prefix for writing files.
    output = output_dir / file_prefix

    # Read options file
    options = Options(output, config)

    if options.parallel_mode:
        parallel_runner(options, output)
    else:
        single_runner(options, output)
