"""
Utility functions for building an empirical quality score model.

The functions in this module extract per-read quality scores from FASTQ files
and compute the statistics required by the ``MarkovQualityModel``.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from ..model_sequencing_error.utils import convert_quality_string
from ..common import open_input

_LOG = logging.getLogger(__name__)

__all__ = [
    "read_quality_lists",
    "compute_initial_distribution",
    "compute_position_distributions",
    "build_markov_model",
]


def read_quality_lists(
    files: Iterable[str],
    max_reads: int,
    offset: int,
) -> Tuple[List[List[int]], int]:
    """Read per-read quality scores from one or more FASTQ files.

    Parameters
    ----------
    files:
        An iterable of FASTQ filenames. Only the first ``max_reads`` reads
        will be consumed from each file. If a file contains fewer than
        ``max_reads`` reads, the entire file is processed.
    max_reads:
        Maximum number of reads per file to consume. A value of ``-1``
        or ``float('inf')`` indicates that all reads in the file should be
        processed.
    offset:
        Numeric quality score offset (usually 33 for Sanger encoding).

    Returns
    -------
    Tuple[List[List[int]], int]
        A tuple containing the list of quality score lists and the read
        length inferred from the first record. If no reads are found this
        function returns an empty list and zero.
    """

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

                # FASTQ format: four lines per record (header, sequence, plus, quality)
                header = fq_in.readline()

                if not header:
                    break  # end of file

                seq = fq_in.readline()
                plus = fq_in.readline()
                qual = fq_in.readline()

                if not qual:
                    break  # truncated file

                quality_scores = convert_quality_string(qual.strip(), offset)

                if not read_length:
                    read_length = len(quality_scores)

                # Skip reads that do not match the inferred read length
                if len(quality_scores) != read_length:
                    _LOG.debug(
                        "Skipping record of length %d (expected %d)",
                        len(quality_scores),
                        read_length,
                    )
                    continue

                qualities.append(quality_scores)
                reads_read += 1

    return qualities, read_length


def compute_initial_distribution(qualities: Iterable[List[int]]) -> Dict[int, float]:
    """Compute the empirical distribution of first quality scores.

    Parameters
    ----------
    qualities:
        An iterable of quality score lists. The first element of each list
        will be used to populate the distribution. If a list is empty, it
        contributes nothing.

    Returns
    -------
    Dict[int, float]
        A dictionary mapping quality scores to counts.
    """

    counts: Dict[int, float] = defaultdict(float)

    for qlist in qualities:
        if not qlist:
            continue

        counts[qlist[0]] += 1.0

    return dict(counts)


def compute_position_distributions(
    qualities: Iterable[List[int]],
    read_length: int,
) -> List[Dict[int, float]]:
    """Compute per-position empirical distributions of quality scores.

    Parameters
    ----------
    qualities:
        An iterable of quality score lists. All lists are expected to have
        length ``read_length``. Any that do not are ignored.
    read_length:
        Length of reads used for training.

    Returns
    -------
    List[Dict[int, float]]
        A list of length ``read_length``. Element ``i`` is a dictionary
        mapping quality scores to counts at position ``i``.
    """

    if read_length <= 0:
        return []

    # One histogram per position
    histograms: List[Dict[int, float]] = [
        defaultdict(float) for _ in range(read_length)
    ]

    for qlist in qualities:

        if len(qlist) != read_length:
            continue

        for i, q in enumerate(qlist):
            histograms[i][int(q)] += 1.0

    # Convert default dicts to plain dicts
    return [dict(h) for h in histograms]


def build_markov_model(
    files: Iterable[str],
    max_reads: int,
    offset: int,
) -> Tuple[Dict[int, float], List[Dict[int, float]], int, int]:
    """Build distributions and max quality for the empirical quality model.

    This is a convenience wrapper around :func:`read_quality_lists`,
    :func:`compute_initial_distribution`, and
    :func:`compute_position_distributions`. It returns the empirical
    distributions and the maximum observed quality, along with the read
    length used for training.

    Parameters
    ----------
    files:
        An iterable of FASTQ filenames.
    max_reads:
        Maximum number of reads per file to consume. ``-1`` or
        ``float('inf')`` indicates that all reads should be processed.
    offset:
        Numeric quality score offset (usually 33 for Sanger encoding).

    Returns
    -------
    Tuple[Dict[int, float], List[Dict[int, float]], int, int]
        A tuple ``(initial_distribution, position_distributions,
        max_quality, read_length)``. Distributions are returned as raw
        counts; the caller is responsible for normalisation.
    """

    qualities, read_length = read_quality_lists(files, max_reads, offset)

    if not qualities:
        raise ValueError("No quality scores could be read from the input files.")

    init_dist = compute_initial_distribution(qualities)
    pos_dists = compute_position_distributions(qualities, read_length)

    # Determine maximum observed quality
    max_quality = 0

    for qlist in qualities:
        if not qlist:
            continue

        mq = max(qlist)

        if mq > max_quality:
            max_quality = mq

    return init_dist, pos_dists, max_quality, read_length
