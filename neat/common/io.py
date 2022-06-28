"""
Functions related to system I/O
"""

__all__ = [
    "is_compressed",
    "open_input",
    "open_output",
    "validate_input_path",
    "validate_output_path"
]

import contextlib
import gzip
import logging
import os
from pathlib import Path
from typing import Callable, Iterator, TextIO

log = logging.getLogger(__name__)


def is_compressed(file: str | Path) -> bool:
    """
    Determine if file is compressed.

    At the moment, function is only able to correctly identify files which were
    gzip compressed

    :param file: Path to a file.

    :return: True if file is compressed, False otherwise.

    Note
    ----
    To determine if the file is gzipped, the function reads the first two bytes of the input file. If these
    bytes are ``1f 8b``, the file is considered to be gzipped as it is highly
    unlikely that an ordinary text files start with those two bytes.
    """
    with open(file, "rb") as buffer:
        magic_number = buffer.read(2)
    if magic_number == b"\x1f\x8b":
        return True
    return False


@contextlib.contextmanager
def open_input(path: str | Path) -> Iterator[TextIO]:
    """
    Opens a file for reading.

    Besides regular text-based files, the function also handles gzipped files.

    :param path: The path to the input file.
    :return: The handle to the text file with input data.
    """
    # Apparently due to a bug mypy doesn't like the ternary operator:
    # - https://github.com/python/mypy/issues/4134
    # - https://stackoverflow.com/a/70832391
    # - https://github.com/python/mypy/issues/12053
    open_: Callable[..., TextIO]
    if is_compressed(path):
        open_ = gzip.open
    else:
        open_ = open
    handle = open_(path, "rt", encoding="utf-8")
    try:
        yield handle
    finally:
        handle.close()


@contextlib.contextmanager
def open_output(path: str | Path, mode: str = 'xt') -> Iterator[TextIO]:
    """
    Opens a file for writing.

    If the directory containing the file does not exist, it will be created
    automatically.

    :param path: The path to the output file.
    :param mode: The mode with which to open the file.
    :return: The handle to the text file where data should be written to.

    Raises
    ------
    FileExistsError
        Raised if the output file already exists.
    PermissionError
        Raised if the calling process does not have adequate access rights to
        write to the output file.
    """
    output_path = Path(path)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    open_: Callable[..., TextIO]
    if 'b' in mode:
        handle = gzip.open(output_path, mode=mode)
    else:
        handle = open(output_path, mode=mode, encoding='utf-8', newline="")

    try:
        yield handle
    finally:
        handle.close()


def validate_input_path(path: str | Path, key: str = None):
    """
    Determine if the input path is valid.

    The input path is valid if it is a file and is not empty. If the path is
    not valid one of the exceptions described later is raised.

    :param path: Path to validate
    :param key: The associated key for the path

    Raises
    ------
    FileNotFoundError
        Raised if the input file does not exist or is not a file.
    RuntimeError
        Raised if the input file is empty.
    PermissionError
        Raised if the calling process has no read access to the file.
    """
    path = Path(path)
    mssg = ''
    if key:
        mssg += f'Invalid value for {key}. '
    if not path.is_file():
        mssg += f"Path '{path}' does not exist or not a file"
        raise FileNotFoundError(mssg)
    stats = path.stat()
    if stats.st_size == 0:
        mssg += f"File '{path}' is empty"
        raise RuntimeError(mssg)
    if not os.access(path, os.R_OK):
        mssg += f"cannot read from '{path}': access denied"
        raise PermissionError(mssg)


def validate_output_path(path: str | Path, is_file: bool = True):
    """
    Determine if the output path is valid.

    If the output file is a directory, it is valid if it exists and the calling
    process has write access to it. If the output path is a file, it is
    valid if the file does not yet exist.

    :param path: The path to validate.
    :param is_file: (optional) If set, validate the path assuming that it points to a file (default).

    Raises
    ------
    FileExistsError
        Raised if path is a file and already exists.
    PermissionError
        Raised if the calling process does not have adequate access rights to.
    """
    path = Path(path)
    if is_file:
        if path.is_file():
            raise FileExistsError(f"file '{path}' already exists")
    else:
        if path.is_dir():
            if not os.access(path, os.W_OK):
                raise PermissionError(f"cannot write to '{path}', access denied")
        else:
            path.mkdir(parents=True)
