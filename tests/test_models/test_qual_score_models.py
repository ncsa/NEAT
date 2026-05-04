"""
Unit tests for MarkovQualityModel
"""

import pytest
import numpy as np
from numpy.random import default_rng

from neat.models.markov_quality_model import MarkovQualityModel

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _simple_model(max_q=42, read_length=151, init_score=30, pos_score=30,
                  transition_distributions=None):
    """Minimal deterministic model with a single quality score at every position."""
    init_dist = {init_score: 1.0}
    pos_dists = [{pos_score: 1.0}] * read_length
    return MarkovQualityModel(init_dist, pos_dists, max_q, read_length,
                              transition_distributions)


# ---------------------------------------------------------------------------
# Original tests — fixed to use valid (non-negative) quality score keys
# ---------------------------------------------------------------------------

def test_markov_quality_model_shapes_and_range():
    """Output array has the requested length and all scores are in [0, max_q]."""
    rng = default_rng(11)
    init_dist = {30: 1.0, 31: 1.0, 32: 1.0}
    pos_dist = {30: 1.0, 31: 2.0, 32: 1.0}   # fix: was {-1, 0, 1} — invalid keys
    max_q = 42
    read_length = 151
    qm = MarkovQualityModel(init_dist, [pos_dist] * read_length, max_q, read_length)
    qs = qm.get_quality_scores(model_read_length=read_length, length=75, rng=rng)
    assert isinstance(qs, np.ndarray)
    assert len(qs) == 75
    assert qs.min() >= 0
    assert qs.max() <= max_q


def test_markov_quality_model_quality_scores_property_matches_range():
    """quality_scores exposes the full contiguous range [0, max_quality]."""
    qm = _simple_model(max_q=40, pos_score=35)
    scores = qm.quality_scores
    assert isinstance(scores, list)
    assert scores[0] == 0
    assert scores[-1] == 40
    assert scores == list(range(0, 41))


def test_markov_quality_model_reproducible_with_seed():
    """Fixed RNG seed produces identical output on two independent calls."""
    init_dist = {30: 1.0, 31: 1.0}
    pos_dist = {30: 1.0, 31: 2.0, 32: 1.0}   # fix: was {-1, 0, 1} — invalid keys
    max_q = 42
    read_length = 151
    qm = MarkovQualityModel(init_dist, [pos_dist] * read_length, max_q, read_length)
    qs1 = qm.get_quality_scores(model_read_length=read_length, length=100, rng=default_rng(12))
    qs2 = qm.get_quality_scores(model_read_length=read_length, length=100, rng=default_rng(12))
    assert np.array_equal(qs1, qs2)


def test_markov_quality_model_invalid_initial_distribution_raises():
    """Empty or zero-mass initial distributions must raise ValueError."""
    read_length = 151
    pos_dists = [{30: 1.0}] * read_length
    with pytest.raises(ValueError):
        MarkovQualityModel({}, pos_dists, 40, read_length)
    with pytest.raises(ValueError):
        MarkovQualityModel({30: 0.0}, pos_dists, 40, read_length)


# ---------------------------------------------------------------------------
# New tests — transition chain path
# ---------------------------------------------------------------------------

def test_markov_model_uses_transition_chain():
    """When transition_distributions is given, scores follow the chain rows."""
    # Transitions always go to 40 regardless of previous score
    always_40 = {q: {40: 1.0} for q in range(0, 43)}
    trans_dists = [always_40] * 150   # read_length - 1
    qm = MarkovQualityModel({30: 1.0}, [{30: 1.0, 40: 1.0}] * 151, 42, 151, trans_dists)
    qs = qm.get_quality_scores(model_read_length=151, length=20, rng=default_rng(0))
    assert qs[0] == 30                              # init distribution
    assert all(qs[i] == 40 for i in range(1, len(qs)))  # chain forces 40


def test_markov_model_falls_back_to_marginal_for_unknown_q_prev():
    """If q_prev has no transition row, the marginal is used without crashing."""
    # Transition map only covers q=99, which is never reached; marginal gives 35
    sparse_trans = [{99: {99: 1.0}}] * 150
    qm = MarkovQualityModel({30: 1.0}, [{35: 1.0}] * 151, 42, 151, sparse_trans)
    qs = qm.get_quality_scores(model_read_length=151, length=10, rng=default_rng(1))
    assert len(qs) == 10
    assert all(0 <= q <= 42 for q in qs)
    # All positions after 0 fall back to marginal → always 35
    assert all(qs[i] == 35 for i in range(1, len(qs)))


def test_markov_model_wrong_transition_length_raises():
    """transition_distributions with wrong length raises ValueError."""
    bad_trans = [{30: {31: 1.0}}] * 50   # should be 150 for read_length=151
    with pytest.raises(ValueError, match="read_length-1"):
        MarkovQualityModel({30: 1.0}, [{30: 1.0}] * 151, 42, 151, bad_trans)


def test_markov_model_output_clipped_to_max_quality():
    """Transition rows that emit out-of-range scores are clipped to max_quality."""
    over_max = {q: {99: 1.0} for q in range(0, 43)}
    trans_dists = [over_max] * 150
    qm = MarkovQualityModel({30: 1.0}, [{99: 1.0}] * 151, 42, 151, trans_dists)
    qs = qm.get_quality_scores(model_read_length=151, length=50, rng=default_rng(5))
    assert qs.max() <= 42


# ---------------------------------------------------------------------------
# New tests — edge cases
# ---------------------------------------------------------------------------

def test_markov_model_length_zero_returns_empty():
    """length=0 must return an empty ndarray without error."""
    qm = _simple_model()
    qs = qm.get_quality_scores(model_read_length=151, length=0, rng=default_rng(0))
    assert isinstance(qs, np.ndarray)
    assert len(qs) == 0


def test_markov_model_length_one_uses_only_init():
    """length=1 reads only from the initial distribution, never marginals."""
    qm = _simple_model(init_score=37, pos_score=10)
    qs = qm.get_quality_scores(model_read_length=151, length=1, rng=default_rng(0))
    assert len(qs) == 1
    assert qs[0] == 37


def test_markov_model_short_read_interpolates_correctly():
    """A read shorter than the model still returns the requested length."""
    qm = _simple_model(read_length=151)
    qs = qm.get_quality_scores(model_read_length=151, length=50, rng=default_rng(0))
    assert len(qs) == 50
    assert all(0 <= q <= 42 for q in qs)


# ---------------------------------------------------------------------------
# New tests — _position_index_for_length
# ---------------------------------------------------------------------------

def test_position_index_for_length_boundaries():
    """First and last positions always map to 0 and read_length-1."""
    qm = _simple_model(read_length=151)
    assert qm._position_index_for_length(0, 50) == 0
    assert qm._position_index_for_length(49, 50) == 150


def test_position_index_for_length_midpoint():
    """Midpoint of a 51-position read maps to midpoint of the 151-position model."""
    qm = _simple_model(read_length=151)
    assert qm._position_index_for_length(25, 51) == 75


def test_position_index_for_length_single_position_read():
    """length=1 always returns index 0 regardless of pos."""
    qm = _simple_model(read_length=151)
    assert qm._position_index_for_length(0, 1) == 0
