"""
Tests for sequencing error model in models
"""
import numpy as np

from Bio.SeqRecord import SeqRecord
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
    reference = SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'), id="fake_1", name="fake", description="fake")
    errors, _padding = model.get_sequencing_errors(padding=5, reference_segment=reference, quality_scores=quality_scores, rng=rng)
    assert all(e.error_type == SingleNucleotideVariant for e in errors)


essentially_zero = 0.0

def test_no_errors_when_avg_zero():
    """When average error is zero, the model should introduce no errors."""
    rng = np.random.default_rng(0)
    model = SequencingErrorModel(read_length=10, avg_seq_error=essentially_zero)
    qs = np.full_like(np.arange(10), 36)
    ref = SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'), id="fake_1", name="fake", description="fake")
    result = model.get_sequencing_errors(padding=5, reference_segment=ref, quality_scores=qs, rng=rng)
    if isinstance(result, tuple):
        errors, pad = result
        assert errors == [] and pad == 5
    else:
        assert result == []