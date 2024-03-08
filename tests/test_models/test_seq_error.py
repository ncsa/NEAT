"""
Tests for sequencing error model in models
"""

import pytest
import numpy as np
import Bio

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from importlib import resources

from neat.models import SequencingErrorModel, default_sequencing_error_model
from neat.variants import SingleNucleotideVariant, Insertion, Deletion


def test_get_seq_error_Insertion():
    """Test that error_type code works correctly"""
    rng = np.random.default_rng(0)
    error_model = SequencingErrorModel(read_length=10, rng=rng, variant_probs={Insertion: 1})
    quality_scores = np.full_like(np.arange(10), 36)
    quality_scores[0] = 0
    reference = SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'), id="fake_1", name="fake", description="fake")
    errors = error_model.get_sequencing_errors(reference, quality_scores)
    assert errors[0][0] == Insertion


def test_get_seq_error_Deletion():
    """Test that error_type code works correctly"""
    rng = np.random.default_rng(0)
    error_model = SequencingErrorModel(read_length=10, rng=rng, variant_probs={Deletion: 1})
    quality_scores = np.full_like(np.arange(10), 36)
    quality_scores[0] = 0
    reference = SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'), id="fake_1", name="fake", description="fake")
    errors = error_model.get_sequencing_errors(reference, quality_scores)
    assert errors[0][0] == Deletion


def test_get_seq_error_SNV():
    """Test that error_type code works correctly"""
    rng = np.random.default_rng(0)
    error_model = SequencingErrorModel(read_length=10, rng=rng, variant_probs={SingleNucleotideVariant: 1})
    quality_scores = np.full_like(np.arange(10), 36)
    quality_scores[0] = 0
    reference = SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'), id="fake_1", name="fake", description="fake")
    errors = error_model.get_sequencing_errors(reference, quality_scores)
    assert errors[0][0] == SingleNucleotideVariant