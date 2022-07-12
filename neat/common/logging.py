"""Sets up logging"""

__all__ = ['setup_logging']

import io
import os
import sys
import logging

from pathlib import Path
from typing import Any

LOG_DETAIL = {
    "LOW": "%(levelname)s:%(name)s:%(message)s",
    "MEDIUM": "%(asctime)s:%(levelname)s:%(name)s:%(message)s",
    "HIGH": "%(asctime)s:%(levelname)s:%(name)s:line %(lineno)s:%(message)s",
}


def setup_logging(
        filename: str,
        severity: str = "INFO",
        verbosity: str = "MEDIUM",
        directory: str = os.getcwd(),
        silent_mode: bool = False
):
    """
    Configure logging for the run

    :param filename: Name to give the log. Default is <timestamp>.neat.log in the current working directory
    :param severity: Severity of events that will be tracked. Defaults to "INFO."
    :param verbosity: Changes the amount of information in the log output.
    :param directory: Directory to store the log. Default is the current working directory.
    :param silent_mode: Default is to output part of the logs to stdout in addition to writing the file
                        Setting this flag will cause it to not print to stdout.
    """

    kwargs: dict[str, Any] = {
        "force": True,
        "format": LOG_DETAIL.get(verbosity.upper(), LOG_DETAIL['MEDIUM'])
    }

    if silent_mode:
        text_trap = io.StringIO()
        sys.stdout = text_trap

    log_file = Path(f'{directory}/{filename}')

    level = getattr(logging, severity.upper(), logging.INFO)
    kwargs['level'] = level

    kwargs['handlers'] = [logging.FileHandler(log_file)]

    logging.basicConfig(**kwargs)


