"""
This function maps the n-gaps for the reference sequence.
"""

import logging
import re

from Bio.Seq import Seq
from typing import Any

__all__ = ['find_n_gaps']

_LOG = logging.getLogger(__name__)


def find_n_gaps(name: str, seq: Seq) -> list[dict[str, Any]]:
    """
    Finds the n-gaps in a given sequence

    :param name: Name of the sequence
    :param seq: Sequence of bases in which to find the N-gaps
    :return: list of dictionaries describing the regions of all N's
    """
    ret_list = []
    for match in re.finditer("N+", str(seq)):
        info = {
            "chrom": name,
            "start": match.start(),
            "end": match.end()
        }
        ret_list.append(info)

    return ret_list
