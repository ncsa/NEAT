"""
Utility functions for building a position-specific Markov quality score model.

The functions in this module extract per-read quality scores from FASTQ files
and compute the information required by the ``MarkovQualityModel``.
"""

import logging
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from ..common import open_input
from ..model_sequencing_error.utils import convert_quality_string

_LOG = logging.getLogger(__name__)

__all__ = [
    "read_quality_lists",
    "compute_initial_distribution",
    "compute_position_distributions",
    "compute_transition_distributions",
    "build_markov_model",
]


def _down_bin_quality(q: int, allowed: List[int]) -> int:
    """
    Map q to the greatest allowed value <= q (down-binning).
    If q is below the smallest allowed, map to allowed[0].
    """

    if not allowed:
        return int(q)

    q = int(q)
    i = bisect_right(allowed, q) - 1

    if i < 0:
        return int(allowed[0])

    return int(allowed[i])


def read_quality_lists(
    files: Iterable[str],
    max_reads: int | float,
    offset: int,
    allowed_quality_scores: Optional[Iterable[int]] = None,
) -> Tuple[List[List[int]], int]:
    """Read per-read quality scores from one or more FASTQ files.

    Returns (qualities, read_length). If no reads are found, returns ([], 0).
    """

    allowed_sorted: Optional[List[int]] = None

    if allowed_quality_scores is not None:
        allowed_sorted = sorted({int(x) for x in allowed_quality_scores})

        if not allowed_sorted:
            allowed_sorted = None

    qualities: List[List[int]] = []
    read_length = 0

    for file in files:

        path = Path(file)
        if not path.exists():
            raise FileNotFoundError(f"Input FASTQ file not found: {file}")

        reads_to_parse = max_reads

        if reads_to_parse in (None, -1):
            reads_to_parse = float("inf")

        reads_read = 0

        with open_input(file) as fq_in:

            while reads_read < reads_to_parse:

                # FASTQ format: four lines per record
                header = fq_in.readline()

                if not header:
                    break  # end of file

                _ = fq_in.readline()  # seq
                _ = fq_in.readline()  # plus
                qual = fq_in.readline()

                if not qual:
                    break

                qlist = convert_quality_string(qual.strip(), offset)

                if not read_length:
                    read_length = len(qlist)

                # Skip reads that do not match the inferred read length
                if len(qlist) != read_length:
                    _LOG.debug(
                        "Skipping record of length %d (expected %d)",
                        len(qlist),
                        read_length,
                    )
                    continue

                if allowed_sorted is not None:
                    qlist = [_down_bin_quality(q, allowed_sorted) for q in qlist]

                qualities.append(qlist)
                reads_read += 1

    return qualities, read_length


def compute_initial_distribution(qualities: Iterable[List[int]]) -> Dict[int, float]:
    """Counts of Q at position 0."""

    counts: Dict[int, float] = defaultdict(float)

    for qlist in qualities:
        if qlist:
            counts[int(qlist[0])] += 1.0

    return dict(counts)


def compute_position_distributions(
    qualities: Iterable[List[int]],
    read_length: int,
) -> List[Dict[int, float]]:
    """Counts of Q at position i."""

    if read_length <= 0:
        return []

    histograms: List[Dict[int, float]] = [defaultdict(float) for _ in range(read_length)]

    for qlist in qualities:

        if len(qlist) != read_length:
            continue

        for i, q in enumerate(qlist):
            histograms[i][int(q)] += 1.0

    return [dict(h) for h in histograms]


def compute_transition_distributions(
        qualities: Iterable[List[int]],
        read_length: int,
) -> List[Dict[int, Dict[int, float]]]:
    """
    Transition counts per position i (0..read_length-2):
      trans[i][q_prev][q_next] += 1
    """

    if read_length <= 1:
        return []

    trans: List[Dict[int, Dict[int, float]]] = [
        defaultdict(lambda: defaultdict(float)) for _ in range(read_length - 1)
    ]

    for qlist in qualities:

        if len(qlist) != read_length:
            continue

        for i in range(read_length - 1):
            q_prev = int(qlist[i])
            q_next = int(qlist[i + 1])
            trans[i][q_prev][q_next] += 1.0

    # Convert nested defaultdicts to plain dicts
    out: List[Dict[int, Dict[int, float]]] = []

    for i in range(read_length - 1):

        pos_dict: Dict[int, Dict[int, float]] = {}

        for q_prev, nexts in trans[i].items():
            pos_dict[int(q_prev)] = {int(qn): float(c) for qn, c in nexts.items()}

        out.append(pos_dict)

    return out


def build_markov_model(
    files: Iterable[str],
    max_reads: int,
    offset: int,
    allowed_quality_scores: Optional[Iterable[int]] = None,
) -> Tuple[
    Dict[int, float],
    List[Dict[int, float]],
    List[Dict[int, Dict[int, float]]],
    int,
    int,
]:
    """Wrapper to create the Markov model."""

    qualities, read_length = read_quality_lists(
        files, max_reads, offset, allowed_quality_scores=allowed_quality_scores
    )

    if not qualities:
        raise ValueError("No quality scores could be read from the input files.")

    init_dist = compute_initial_distribution(qualities)
    pos_dists = compute_position_distributions(qualities, read_length)
    trans_dists = compute_transition_distributions(qualities, read_length)

    # Determine maximum quality (post-binning, if applied)
    max_quality = 0
    for qlist in qualities:
        if qlist:
            max_quality = max(max_quality, max(qlist))

    # If user supplied bins, max_quality should not exceed the max bin.
    if allowed_quality_scores is not None:
        bins = sorted({int(x) for x in allowed_quality_scores})

        if bins:
            max_quality = min(max_quality, max(bins))

    return init_dist, pos_dists, trans_dists, int(max_quality), int(read_length)
