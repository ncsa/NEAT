"""
Unit tests for MarkovQualityModel
"""

import pytest
import numpy as np
from numpy.random import default_rng

from neat.models.markov_quality_model import MarkovQualityModel


def test_markov_quality_model_shapes_and_range():
    """
    Basic sanity check for MarkovQualityModel: output array shape and bounds.
    """
    rng = default_rng(11)
    # Simple symmetric distributions around a high-quality region
    init_dist = {30: 1.0, 31: 1.0, 32: 1.0}
    single_pos_dist = {-1: 1.0, 0: 2.0, 1: 1.0}
    max_q = 42
    read_length = 151
    # Position-specific distributions that reuse the same shape at every position
    position_distributions = [single_pos_dist] * read_length
    qm = MarkovQualityModel(init_dist, position_distributions, max_q, read_length)
    qs = qm.get_quality_scores(model_read_length=read_length, length=75, rng=rng)
    assert isinstance(qs, np.ndarray)
    assert len(qs) == 75
    # Ensure we stay within the valid range
    assert qs.min() >= 0
    assert qs.max() <= max_q


def test_markov_quality_model_quality_scores_property_matches_range():
    """
    quality_scores should expose the full discrete range [0, max_quality].
    """
    init_dist = {35: 1.0}
    single_pos_dist = {0: 1.0}
    max_q = 40
    read_length = 151
    position_distributions = [single_pos_dist] * read_length
    qm = MarkovQualityModel(init_dist, position_distributions, max_q, read_length)
    scores = qm.quality_scores
    assert isinstance(scores, list)
    assert scores[0] == 0
    assert scores[-1] == max_q
    # The range should be contiguous
    assert scores == list(range(0, max_q + 1))


def test_markov_quality_model_reproducible_with_seed():
    """Markov quality model should be deterministic for a fixed RNG state."""
    init_dist = {30: 1.0, 31: 1.0}
    single_pos_dist = {-1: 1.0, 0: 2.0, 1: 1.0}
    max_q = 42
    read_length = 151
    position_distributions = [single_pos_dist] * read_length
    qm = MarkovQualityModel(init_dist, position_distributions, max_q, read_length)
    rng1 = default_rng(12)
    rng2 = default_rng(12)
    qs1 = qm.get_quality_scores(model_read_length=read_length, length=100, rng=rng1)
    qs2 = qm.get_quality_scores(model_read_length=read_length, length=100, rng=rng2)
    assert isinstance(qs1, np.ndarray)
    assert isinstance(qs2, np.ndarray)
    assert np.array_equal(qs1, qs2)


def test_markov_quality_model_invalid_initial_distribution_raises():
    """
    The model should reject empty or zero-mass initial distributions.
    """
    read_length = 151
    single_pos_dist = {0: 1.0}
    position_distributions = [single_pos_dist] * read_length
    # Empty initial distribution
    with pytest.raises(ValueError):
        MarkovQualityModel({}, position_distributions, 40, read_length)
    # Zero total mass in initial distribution
    with pytest.raises(ValueError):
        MarkovQualityModel({30: 0.0}, position_distributions, 40, read_length)
