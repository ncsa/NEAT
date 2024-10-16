"""
Creates a model of quality scores by position along the read.
"""

import gzip
import pickle
import numpy as np
import logging

from pathlib import Path

from .utils import make_qual_score_list, apply_markov_chain, save_file
from ..models import QualityScoreModel

__all__ = [
    "quality_score_model_runner"
]


def quality_score_model_runner(
        input_file: str,
        output_path: str,
):
    """
    Sets up and calls the sequencing error modeling core code.

    :param input_file: This is a fastq input file for quality score modeling (preferred over a bam file)
    :param output_path: This is the primary output that stores the Markov model's predictions.
    """

    quality_df = make_qual_score_list(input_file)
    markov_preds_df = apply_markov_chain(quality_df)
    save_file(markov_preds_df, output_path, output_path)
