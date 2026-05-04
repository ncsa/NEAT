"""
Tests for sequencing error model in models
"""

import numpy as np

from Bio.Seq import Seq

from neat.models import SequencingErrorModel
from neat.variants import SingleNucleotideVariant, Insertion, Deletion

def test_get_seq_error_snv_only():
    """With SNV-only probabilities, produced errors should be SNVs."""
    rng = np.random.default_rng(0)
    model = SequencingErrorModel(read_length=10, variant_probs={SingleNucleotideVariant: 1.0})
    quality_scores = np.full_like(np.arange(10), 36)
    # Make at least one base very error-prone
    quality_scores[0] = 0
    reference = Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    errors, _padding = model.get_sequencing_errors(padding=5, reference_segment=reference, quality_scores=quality_scores, num_errors=3, rng=rng)
    assert all(e.error_type == SingleNucleotideVariant for e in errors)


    # test_no_errors_when_avg_zero removed:
    # duplicate of test_error_models.py::test_sem_zero_error_rate_returns_empty