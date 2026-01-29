"""
Runner for creating quality score models.

This module implements the core logic for generating quality score models
from input FASTQ files.

The resulting models are saved as a gzip-compressed pickle file.
"""

import gzip
import pickle
import logging
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np

from ..common import validate_input_path, validate_output_path
from ..models import SequencingErrorModel, TraditionalQualityModel
from ..model_sequencing_error.utils import parse_file
from ..models.markov_quality_model import MarkovQualityModel
from ..quality_score_modeling.markov_utils import build_markov_model

__all__ = ["model_qual_score_runner"]

_LOG = logging.getLogger(__name__)


def _prepare_quality_scores_argument(
    qual_scores: Iterable[int | float],
) -> Tuple[List[int], bool]:
    """Normalize the quality scores argument to a list of ints."""

    if isinstance(qual_scores, int):
        max_q = int(qual_scores)
        return list(range(0, max_q + 1)), False

    sorted_scores = sorted(int(x) for x in qual_scores)
    max_q = max(sorted_scores)

    return list(range(0, max_q + 1)), True


def model_qual_score_runner(
    files: List[str],
    offset: int,
    qual_scores: Iterable[int | float],
    max_reads: int,
    use_markov: bool,
    overwrite: bool,
    output_dir: str,
    output_prefix: str,
) -> None:
    """Create and save a quality score model from FASTQ data.

    Parameters
    ----------
    files:
        A list of FASTQ file names.
    offset:
        Quality score offset (e.g. 33 for Sanger encoding).
    qual_scores:
        Either a single integer specifying the maximum quality score or an
        iterable of integers specifying the allowed quality scores (binned
        model).
    max_reads:
        Maximum number of reads to process from each input file. A value of
        ``-1`` or ``float('inf')`` processes all reads.
    use_markov:
        If ``True``, construct a Markov chain–based quality model. If
        ``False``, only the traditional quality model is produced.
    overwrite:
        If ``True``, allow an existing output file to be overwritten.
    output_dir:
        Directory to write the output model file.
    output_prefix:
        Prefix to use for the output file name. The file will be named
        ``{output_prefix}.p.gz`` in the ``output_dir``.
    """

    if len(files) > 2:
        _LOG.info("Only processing the first two input files")
        files = files[:2]

    # Validate input paths
    for file in files:
        validate_input_path(file)

    _LOG.debug("Input files: %s", ", ".join(str(x) for x in files))
    _LOG.debug("Quality offset: %d", offset)

    final_quality_scores, binned = _prepare_quality_scores_argument(qual_scores)
    _LOG.debug("Quality scores: %s", final_quality_scores)

    # Determine maximum number of reads to process
    if max_reads in (-1, None):
        num_records_to_process = float("inf")

    else:
        num_records_to_process = max_reads
    _LOG.debug(
        "Maximum number of records to process: %s",
        "all" if num_records_to_process == float("inf") else num_records_to_process,
    )

    # Validate output directory and file
    validate_output_path(output_dir, is_file=False)
    output_path = Path(output_dir)
    output_file = output_path / f"{output_prefix}.p.gz"
    validate_output_path(output_file, overwrite=overwrite)
    _LOG.info("Writing output to: %s", output_file)

    # Containers for per-file quality model parameters
    read_parameters: List[np.ndarray] = []
    average_errors: List[float] = []
    read_length = 0

    # Collect traditional model parameters using existing NEAT utilities
    file_num = 0

    for file in files:

        file_num += 1
        _LOG.info("Reading file %d of %d", file_num, len(files))

        params_by_position, file_avg_error, read_length = parse_file(
            file,
            final_quality_scores,
            num_records_to_process,
            offset,
            read_length,
        )

        read_parameters.append(params_by_position)
        average_errors.append(file_avg_error)
        _LOG.info("Finished reading file %d", file_num)

    if not read_parameters:
        raise RuntimeError(
            "No quality score parameters were computed; check input FASTQ files."
        )

    average_error = float(np.average(average_errors)) if average_errors else 0.0
    _LOG.info("Average sequencing error across files: %f", average_error)

    # Prepare models for each input file
    models: List[Tuple[SequencingErrorModel, TraditionalQualityModel, Optional[MarkovQualityModel]]] = []

    for idx in range(len(read_parameters)):

        # Sequencing error model (always produced)
        seq_err_model = SequencingErrorModel(
            avg_seq_error=average_error,
            read_length=read_length,
        )

        # Traditional quality model (always produced)
        trad_model = TraditionalQualityModel(
            average_error=average_error,
            quality_scores=np.array(final_quality_scores),
            qual_score_probs=read_parameters[idx],
        )

        # Optionally build Markov quality model
        markov_model: Optional[MarkovQualityModel] = None

        if use_markov:
            try:
                (
                    init_dist,
                    pos_dists,
                    max_quality,
                    train_read_length,
                ) = build_markov_model(
                    [files[idx]],
                    num_records_to_process,
                    offset,
                )

                markov_model = MarkovQualityModel(
                    initial_distribution=init_dist,
                    position_distributions=pos_dists,
                    max_quality=max_quality,
                    read_length=train_read_length,
                )

            except Exception as exc:
                _LOG.error(
                    "Failed to construct Markov model for %s: %s",
                    files[idx],
                    exc,
                )
                raise

        models.append((seq_err_model, trad_model, markov_model))

    _LOG.debug("Constructed %d model(s)", len(models))

    # Write out the models
    with gzip.open(output_file, "wb") as out_model:
        if not models:
            raise RuntimeError(
                "Internal error: no models constructed, so there is nothing to write."
            )

        if len(models) == 1:
            seq_err1, trad1, markov1 = models[0]
            pickle.dump(
                {
                    "error_model1": seq_err1,
                    "error_model2": None,
                    "qual_score_model1": markov1 if use_markov else trad1,
                    "qual_score_model2": None,
                },
                out_model,
            )

        elif len(models) == 2:
            (seq_err1, trad1, markov1), (seq_err2, trad2, markov2) = models
            pickle.dump(
                {
                    "error_model1": seq_err1,
                    "error_model2": seq_err2,
                    "qual_score_model1": markov1 if use_markov else trad1,
                    "qual_score_model2": markov2 if use_markov else trad2,
                },
                out_model,
            )

        else:
            # NEAT's read simulator only understands one or two models
            raise RuntimeError(
                f"Expected at most two quality models, but constructed {len(models)}."
            )

    _LOG.info("Quality score model saved to %s", output_file)
