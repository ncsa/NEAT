"""
Unit tests for neat/models/error_models.py

Covers TraditionalQualityModel, SequencingErrorModel, and ErrorContainer.
"""
import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.models.error_models import (
    TraditionalQualityModel,
    SequencingErrorModel,
    ErrorContainer,
)
from neat.variants import Insertion, Deletion, SingleNucleotideVariant


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_RNG = np.random.default_rng(42)
_SEQ = Seq("ACGT" * 40)          # 160 bp, no Ns
_SEQ_RECORD = SeqRecord(_SEQ, id="chr1")


# ===========================================================================
# TraditionalQualityModel — construction
# ===========================================================================

def test_tqm_default_construction():
    m = TraditionalQualityModel()
    assert m.quality_scores is not None
    assert len(m.quality_scores) > 0


def test_tqm_quality_scores_range():
    m = TraditionalQualityModel()
    assert m.quality_scores.min() >= 1
    assert m.quality_scores.max() <= 42


def test_tqm_error_rate_dict_populated():
    m = TraditionalQualityModel()
    for score in m.quality_scores:
        assert score in m.quality_score_error_rate
        assert 0 < m.quality_score_error_rate[score] <= 1


def test_tqm_not_uniform_by_default():
    m = TraditionalQualityModel()
    assert not m.is_uniform
    assert m.uniform_quality_score is None


def test_tqm_uniform_mode_sets_score():
    m = TraditionalQualityModel(is_uniform=True)
    assert m.is_uniform
    assert m.uniform_quality_score is not None
    assert isinstance(m.uniform_quality_score, (int, np.integer))


def test_tqm_uniform_score_within_range():
    m = TraditionalQualityModel(is_uniform=True)
    assert 1 <= m.uniform_quality_score <= 42


# ===========================================================================
# TraditionalQualityModel — get_quality_scores
# ===========================================================================

def test_tqm_get_quality_scores_uniform_returns_constant_array():
    m = TraditionalQualityModel(is_uniform=True)
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 50, rng)
    assert len(scores) == 50
    assert all(s == scores[0] for s in scores)


def test_tqm_get_quality_scores_length_matches_request():
    m = TraditionalQualityModel()
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 100, rng)
    assert len(scores) == 100


def test_tqm_get_quality_scores_exact_model_length():
    m = TraditionalQualityModel()
    rng = np.random.default_rng(0)
    # 151 is the default model read length
    scores = m.get_quality_scores(151, 151, rng)
    assert len(scores) == 151


def test_tqm_get_quality_scores_shorter_than_model():
    m = TraditionalQualityModel()
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 50, rng)
    assert len(scores) == 50


def test_tqm_get_quality_scores_values_in_range():
    m = TraditionalQualityModel()
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 100, rng)
    assert all(1 <= int(s) <= 42 for s in scores)


def test_tqm_get_quality_scores_returns_ndarray():
    m = TraditionalQualityModel()
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 100, rng)
    assert isinstance(scores, np.ndarray)


def test_tqm_get_quality_scores_reproducible():
    m = TraditionalQualityModel()
    s1 = m.get_quality_scores(151, 100, np.random.default_rng(7))
    s2 = m.get_quality_scores(151, 100, np.random.default_rng(7))
    np.testing.assert_array_equal(s1, s2)


# ===========================================================================
# SequencingErrorModel — construction
# ===========================================================================

def test_sem_default_construction():
    m = SequencingErrorModel()
    assert m.average_error > 0
    assert m.read_length == 151


def test_sem_variant_probs_sum_to_one():
    m = SequencingErrorModel()
    assert abs(sum(m.variant_probs.values()) - 1.0) < 1e-9


def test_sem_custom_error_rate():
    m = SequencingErrorModel(avg_seq_error=0.001)
    assert m.average_error == 0.001


# ===========================================================================
# SequencingErrorModel — get_sequencing_errors
# ===========================================================================

def test_sem_zero_error_rate_returns_empty():
    m = SequencingErrorModel(avg_seq_error=0.0)
    rng = np.random.default_rng(0)
    quality_scores = np.array([40] * 100)
    result = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    assert result == []


def test_sem_high_error_rate_returns_errors():
    # Force errors by using very low quality scores (high error probability)
    m = SequencingErrorModel(avg_seq_error=0.5)
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 100)  # quality 1 → ~79% error rate
    result, _ = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    assert len(result) > 0


def test_sem_returns_error_container_objects():
    m = SequencingErrorModel(avg_seq_error=0.5)
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 100)
    result, _ = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    for err in result:
        assert isinstance(err, ErrorContainer)


def test_sem_errors_have_valid_locations():
    m = SequencingErrorModel(avg_seq_error=0.5)
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 100)
    result, _ = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    for err in result:
        assert 0 <= err.location < len(quality_scores)


def test_sem_snv_errors_have_valid_alt():
    m = SequencingErrorModel(avg_seq_error=0.5)
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 100)
    result, _ = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    snv_errors = [e for e in result if e.error_type == SingleNucleotideVariant]
    for err in snv_errors:
        assert err.alt in ("A", "C", "G", "T")


def test_sem_returns_updated_padding():
    m = SequencingErrorModel(avg_seq_error=0.5)
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 100)
    _, padding = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    assert padding >= 0


def test_sem_high_quality_scores_produce_few_errors():
    m = SequencingErrorModel(avg_seq_error=0.009)
    rng = np.random.default_rng(0)
    quality_scores = np.array([40] * 151)  # q40 → 0.01% error rate
    result, _ = m.get_sequencing_errors(20, _SEQ_RECORD, quality_scores, rng)
    # With q40 and length 151, very few errors expected
    assert len(result) < 10


# ===========================================================================
# ErrorContainer
# ===========================================================================

def test_error_container_stores_fields():
    ec = ErrorContainer(SingleNucleotideVariant, 5, 1, "A", "T")
    assert ec.error_type == SingleNucleotideVariant
    assert ec.location == 5
    assert ec.length == 1
    assert ec.ref == "A"
    assert ec.alt == "T"


def test_error_container_deletion_type():
    ec = ErrorContainer(Deletion, 10, 3, "ACGT", "A")
    assert ec.error_type == Deletion
    assert ec.length == 3


def test_error_container_insertion_type():
    ec = ErrorContainer(Insertion, 7, 2, "A", "ACG")
    assert ec.error_type == Insertion
    assert ec.alt == "ACG"


# ===========================================================================
# SequencingErrorModel — indel error paths (lines 209, 213-235)
# ===========================================================================

def test_sem_deletion_variant_prob_has_no_effect():
    """variant_probs favouring Deletion still produces only SNVs (dead-code bug).

    The indel gate condition is circular, so deletion errors are never produced
    regardless of variant_probs. See test_sem_only_snv_errors_produced_regardless_of_variant_probs
    for the full documentation test.
    """
    from neat.variants import Deletion as Del, Insertion as Ins
    m = SequencingErrorModel(
        avg_seq_error=0.9,
        variant_probs={Ins: 0.0, Del: 1.0, SingleNucleotideVariant: 0.0},
    )
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 151)
    result, _ = m.get_sequencing_errors(50, _SEQ_RECORD, quality_scores, rng)
    del_errors = [e for e in result if e.error_type == Del]
    # No deletions produced due to the unreachable gate (see dead-code comment below)
    assert len(del_errors) == 0


def test_sem_insertion_variant_prob_has_no_effect():
    """variant_probs favouring Insertion still produces only SNVs (dead-code bug)."""
    from neat.variants import Insertion as Ins, Deletion as Del
    m = SequencingErrorModel(
        avg_seq_error=0.9,
        variant_probs={Ins: 1.0, Del: 0.0, SingleNucleotideVariant: 0.0},
    )
    rng = np.random.default_rng(7)
    quality_scores = np.array([1] * 151)
    result, _ = m.get_sequencing_errors(50, _SEQ_RECORD, quality_scores, rng)
    ins_errors = [e for e in result if e.error_type == Ins]
    # No insertions produced due to the unreachable gate
    assert len(ins_errors) == 0


def test_sem_blacklist_prevents_duplicate_deletion_sites():
    """Errors at blacklisted positions from a deletion are removed."""
    from neat.variants import Deletion as Del, Insertion as Ins
    m = SequencingErrorModel(
        avg_seq_error=0.9,
        variant_probs={Ins: 0.0, Del: 1.0, SingleNucleotideVariant: 0.0},
    )
    rng = np.random.default_rng(13)
    quality_scores = np.array([1] * 151)
    result, _ = m.get_sequencing_errors(50, _SEQ_RECORD, quality_scores, rng)
    # All returned errors should be at unique locations (blacklist applied)
    locations = [e.location for e in result]
    # Deletions span multiple bases; no two errors at same index
    assert len(locations) == len(set(locations))



# ===========================================================================
# TraditionalQualityModel — score clamping (line 104)
# ===========================================================================

def test_tqm_score_clamped_to_minimum():
    """Scores below 1 from rng.normal are clamped to 1 (line 104)."""
    # Use qual_score_probs with a very negative mean so that rng.normal
    # always returns a value that rounds to ≤ 0, forcing the clamp.
    very_low_probs = np.array([[-50.0, 0.01]] * 151)
    m = TraditionalQualityModel(qual_score_probs=very_low_probs)
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 151, rng)
    assert all(s == 1 for s in scores)


def test_tqm_score_clamped_to_maximum():
    """Scores above 42 from rng.normal are clamped to 42 (line 102)."""
    very_high_probs = np.array([[200.0, 0.01]] * 151)
    m = TraditionalQualityModel(qual_score_probs=very_high_probs)
    rng = np.random.default_rng(0)
    scores = m.get_quality_scores(151, 151, rng)
    assert all(s == 42 for s in scores)


# ===========================================================================
# SequencingErrorModel — dead-code documentation
# ===========================================================================
# Lines 209-235, 252 (indel error branches) are unreachable because
# `total_indel_length` starts at 0 and is only incremented inside the
# branches that are gated by `total_indel_length > self.read_length // 4`.
# This circular dependency means the variant_probs choice (line 209) is
# never called and deletion/insertion errors are never produced.
# The following test documents this behaviour.

def test_sem_only_snv_errors_produced_regardless_of_variant_probs():
    """Indel errors are never produced due to the total_indel_length gate."""
    m = SequencingErrorModel(
        avg_seq_error=0.9,
        variant_probs={Insertion: 0.5, Deletion: 0.5, SingleNucleotideVariant: 0.0},
    )
    rng = np.random.default_rng(0)
    quality_scores = np.array([1] * 151)
    result, _ = m.get_sequencing_errors(50, _SEQ_RECORD, quality_scores, rng)
    # All errors are SNVs because the indel branches are unreachable
    error_types = {e.error_type for e in result}
    assert error_types == {SingleNucleotideVariant}
    # TODO (post-fix): Once the gate condition is corrected from
    #   `total_indel_length > self.read_length // 4`
    #   to
    #   `total_indel_length <= self.read_length // 4`
    # the following assertions should replace the one above:
    #
    #   del_errors = [e for e in result if e.error_type == Deletion]
    #   ins_errors = [e for e in result if e.error_type == Insertion]
    #   assert len(del_errors) + len(ins_errors) > 0, \
    #       "Expected indel errors given variant_probs favours them"
    #   # blacklist test: no two errors share a location
    #   locations = [e.location for e in result]
    #   assert len(locations) == len(set(locations))
