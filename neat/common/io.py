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
from Bio import bgzf

log = logging.getLogger(__name__)


def is_compressed(file: str | Path) -> bool:
    """
    Determine if file_list is compressed.

    At the moment, function is only able to correctly identify files which were
    gzip compressed

    :param file: Path to a file_list.

    :return: True if file_list is compressed, False otherwise.

    Note
    ----
    To determine if the file_list is gzipped, the function reads the first two bytes of the input file_list. If these
    bytes are ``1f 8b``, the file_list is considered to be gzipped as it is highly
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
    Opens a file_list for reading.

    Besides regular text-based files, the function also handles gzipped files.

    :param path: The path to the input file_list.
    :return: The handle to the text file_list with input data.
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
def open_output(path: str | Path, mode: str = 'wt') -> Iterator[TextIO]:
    """
    Opens a file_list for writing.

    If the directory containing the file_list does not exist, it will be created
    automatically.

    :param path: The path to the output file_list.
    :param mode: The mode with which to open the file_list.
    :return: The handle to the text file_list where data should be written to.

    Raises
    ------
    FileExistsError
        Raised if the output file_list already exists.
    PermissionError
        Raised if the calling process does not have adequate access rights to
        write to the output file_list.
    """
    output_path = Path(path)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    open_: Callable[..., TextIO]
    # This is a set intersection to see if the file has gz or bgz in name. This list can be expanded as needed.
    if {'.gz', '.bgz'} & set(output_path.suffixes):
        open_ = bgzf.open
    else:
        open_ = open
    handle = open_(output_path, mode=mode)

    try:
        yield handle
    finally:
        handle.close()


def validate_input_path(path: str | Path, key: str = None):
    """
    Determine if the input path is valid.

    The input path is valid if it is a file_list and is not empty. If the path is
    not valid one of the exceptions described later is raised.

    :param path: Path to validate
    :param key: The associated key for the path

    Raises
    ------
    FileNotFoundError
        Raised if the input file_list does not exist or is not a file_list.
    RuntimeError
        Raised if the input file_list is empty.
    PermissionError
        Raised if the calling process has no read access to the file_list.
    """
    path = Path(path)
    mssg = ''
    if key:
        mssg += f'Invalid value for {key}. '
    if not path.is_file():
        mssg += f"Path '{path}' does not exist or not a file_list"
        raise FileNotFoundError(mssg)
    stats = path.stat()
    if stats.st_size == 0:
        mssg += f"File '{path}' is empty"
        raise RuntimeError(mssg)
    if not os.access(path, os.R_OK):
        mssg += f"cannot read from '{path}': access denied"
        raise PermissionError(mssg)


def validate_output_path(path: str | Path, is_file: bool = True, overwrite: bool = False):
    """
    Determine if the output path is valid.

    If the output file_list is a directory, it is valid if it exists and the calling
    process has write access to it. If the output path is a file_list, it is
    valid if the file_list does not yet exist.

    :param path: The path to validate.
    :param is_file: (optional) If set, validate the path assuming that it points to a file_list (default).
    :param overwrite: (optional) If set, NEAT will overwrite existing output

    Raises
    ------
    FileExistsError
        Raised if path is a file_list and already exists.
    PermissionError
        Raised if the calling process does not have adequate access rights to.
    """
    path = Path(path)
    if is_file:
        if path.is_file() and not overwrite:
            raise FileExistsError(f"file_list '{path}' already exists")
    else:
        if path.is_dir():
            if not os.access(path, os.W_OK):
                raise PermissionError(f"cannot write to '{path}', access denied")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
