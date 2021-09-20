"""
This file contains several standard functions that will be used throughout the program. Each function checks input
and issues an error if there is something wrong.
"""

import pathlib
import sys
import logging
from source.error_handling import will_exit


def check_file_open(filename: str, err_string: str, required: bool = False) -> None:
    """
    Checks that the filename is not empty and that it is indeed a  file

    :param filename: file name, string
    :param err_string: string of the error if it is not a file
    :param required: If not required, skips the check
    :return: None
    """
    if required or filename is not None:
        if filename is None:
            print('\n' + err_string + '\n')
            logging.error(err_string)
            will_exit(1)
        else:
            try:
                pathlib.Path(filename).resolve(strict=True)
            except FileNotFoundError:
                print('\n' + err_string + '\n')
                logging.error(err_string)
                will_exit(1)


def is_in_range(value: float, lower_bound: float, upper_bound: float, err_string: str) -> None:
    """
    Checks that value is between the lower bound and upper bound, and if not prints an error message
    (err_string) and exits the program.

    :param value: float for the value
    :param lower_bound: float for the upper bound
    :param upper_bound: float for the lower bound
    :param err_string: string of the error message to print if the value is out of range
    :return: None
    """
    if value < lower_bound or value > upper_bound:
        print('\n' + err_string + '\n')
        logging.error(err_string)
        will_exit(1)
