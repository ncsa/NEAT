"""
Unit tests for neat/quality_score_modeling/markov_utils.py
"""

import pytest

from neat.quality_score_modeling.markov_utils import (
    _down_bin_quality,
    read_quality_lists,
    compute_initial_distribution,
    compute_position_distributions,
    compute_transition_distributions,
    build_markov_model,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fastq(path, reads):
    """Write a list of (seq, qual_string) tuples as a FASTQ file."""
    lines = []
    for i, (seq, qual) in enumerate(reads):
        lines += [f"@read{i}", seq, "+", qual]
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# _down_bin_quality
# ---------------------------------------------------------------------------

def test_down_bin_exact_match():
    assert _down_bin_quality(30, [10, 20, 30, 40]) == 30


def test_down_bin_between_bins_maps_down():
    assert _down_bin_quality(25, [10, 20, 30, 40]) == 20


def test_down_bin_below_min_maps_to_first_bin():
    assert _down_bin_quality(5, [10, 20, 30]) == 10


def test_down_bin_above_max_maps_to_last_bin():
    assert _down_bin_quality(99, [10, 20, 30]) == 30


def test_down_bin_empty_allowed_returns_q_unchanged():
    assert _down_bin_quality(25, []) == 25


def test_down_bin_single_bin():
    assert _down_bin_quality(0, [20]) == 20
    assert _down_bin_quality(20, [20]) == 20
    assert _down_bin_quality(40, [20]) == 20


# ---------------------------------------------------------------------------
# compute_initial_distribution
# ---------------------------------------------------------------------------

def test_compute_initial_distribution_basic():
    quals = [[30, 31, 32], [30, 28, 27], [35, 30, 29]]
    result = compute_initial_distribution(quals)
    assert result[30] == 2.0
    assert result[35] == 1.0
    assert 28 not in result   # position 0 only
    assert 31 not in result


def test_compute_initial_distribution_empty_input():
    assert compute_initial_distribution([]) == {}


def test_compute_initial_distribution_skips_empty_reads():
    result = compute_initial_distribution([[], [30, 31]])
    assert result == {30: 1.0}


def test_compute_initial_distribution_single_read():
    result = compute_initial_distribution([[40, 38, 35]])
    assert result == {40: 1.0}


# ---------------------------------------------------------------------------
# compute_position_distributions
# ---------------------------------------------------------------------------

def test_compute_position_distributions_basic():
    quals = [[30, 35, 40], [30, 36, 41]]
    result = compute_position_distributions(quals, 3)
    assert len(result) == 3
    assert result[0][30] == 2.0
    assert result[1][35] == 1.0
    assert result[1][36] == 1.0
    assert result[2][40] == 1.0
    assert result[2][41] == 1.0


def test_compute_position_distributions_skips_length_mismatch():
    quals = [[30, 31, 32], [30, 31]]   # second read is wrong length
    result = compute_position_distributions(quals, 3)
    assert result[0][30] == 1.0   # only first read counted
    assert 31 not in result[0]


def test_compute_position_distributions_zero_length_returns_empty():
    assert compute_position_distributions([[30, 31]], 0) == []


def test_compute_position_distributions_single_position():
    quals = [[40], [38], [40]]
    result = compute_position_distributions(quals, 1)
    assert len(result) == 1
    assert result[0][40] == 2.0
    assert result[0][38] == 1.0


# ---------------------------------------------------------------------------
# compute_transition_distributions
# ---------------------------------------------------------------------------

def test_compute_transition_distributions_basic():
    quals = [[30, 31, 32], [30, 31, 33]]
    result = compute_transition_distributions(quals, 3)
    assert len(result) == 2                  # read_length - 1
    assert result[0][30][31] == 2.0          # 30→31 seen twice at position 0
    assert result[1][31][32] == 1.0          # 31→32 once at position 1
    assert result[1][31][33] == 1.0          # 31→33 once at position 1


def test_compute_transition_distributions_read_length_one_returns_empty():
    assert compute_transition_distributions([[30]], 1) == []


def test_compute_transition_distributions_read_length_zero_returns_empty():
    assert compute_transition_distributions([], 0) == []


def test_compute_transition_distributions_skips_length_mismatch():
    quals = [[30, 31, 32], [30, 31]]   # second read is wrong length
    result = compute_transition_distributions(quals, 3)
    assert result[0][30][31] == 1.0    # only first read counted


def test_compute_transition_distributions_self_transitions():
    """A read with a constant quality score produces only self-transitions."""
    quals = [[35, 35, 35, 35]]
    result = compute_transition_distributions(quals, 4)
    assert len(result) == 3
    for pos in result:
        assert pos[35][35] == 1.0
        assert len(pos) == 1


# ---------------------------------------------------------------------------
# read_quality_lists
# ---------------------------------------------------------------------------

def test_read_quality_lists_basic(tmp_path):
    fq = tmp_path / "test.fastq"
    # 'I' = ASCII 73, offset 33 → score 40
    _write_fastq(fq, [("ACGTA", "IIIII"), ("ACGTA", "IIIII")])
    quals, read_length = read_quality_lists([str(fq)], max_reads=100, offset=33)
    assert read_length == 5
    assert len(quals) == 2
    assert quals[0] == [40, 40, 40, 40, 40]


def test_read_quality_lists_applies_offset(tmp_path):
    fq = tmp_path / "test.fastq"
    # '!' = ASCII 33, offset 33 → score 0
    _write_fastq(fq, [("ACGT", "!!!!")])
    quals, _ = read_quality_lists([str(fq)], max_reads=100, offset=33)
    assert quals[0] == [0, 0, 0, 0]


def test_read_quality_lists_respects_max_reads(tmp_path):
    fq = tmp_path / "test.fastq"
    _write_fastq(fq, [("ACGT", "IIII")] * 10)
    quals, _ = read_quality_lists([str(fq)], max_reads=3, offset=33)
    assert len(quals) == 3


def test_read_quality_lists_applies_binning(tmp_path):
    fq = tmp_path / "test.fastq"
    # Scores 40 ('I'), binned to allowed [10, 20, 30] → maps to 30
    _write_fastq(fq, [("ACGT", "IIII")])
    quals, _ = read_quality_lists([str(fq)], max_reads=100, offset=33,
                                  allowed_quality_scores=[10, 20, 30])
    assert quals[0] == [30, 30, 30, 30]


def test_read_quality_lists_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        read_quality_lists(["/nonexistent/path.fastq"], max_reads=10, offset=33)


def test_read_quality_lists_empty_file_returns_empty(tmp_path):
    fq = tmp_path / "empty.fastq"
    fq.write_text("")
    quals, read_length = read_quality_lists([str(fq)], max_reads=100, offset=33)
    assert quals == []
    assert read_length == 0


# ---------------------------------------------------------------------------
# build_markov_model (end-to-end)
# ---------------------------------------------------------------------------

def test_build_markov_model_basic(tmp_path):
    fq = tmp_path / "test.fastq"
    _write_fastq(fq, [("ACGTA", "IIIII"), ("ACGTA", "IIIII")])
    init, pos_dists, trans_dists, max_q, read_len = build_markov_model(
        [str(fq)], max_reads=100, offset=33
    )
    assert read_len == 5
    assert max_q == 40
    assert init == {40: 2.0}
    assert len(pos_dists) == 5
    assert len(trans_dists) == 4   # read_length - 1
    assert trans_dists[0][40][40] == 2.0


def test_build_markov_model_max_quality_capped_by_bins(tmp_path):
    fq = tmp_path / "test.fastq"
    _write_fastq(fq, [("ACGT", "IIII")])   # score 40
    _, _, _, max_q, _ = build_markov_model(
        [str(fq)], max_reads=100, offset=33,
        allowed_quality_scores=[10, 20, 30]
    )
    assert max_q == 30   # capped to max allowed bin


def test_build_markov_model_no_reads_raises(tmp_path):
    fq = tmp_path / "empty.fastq"
    fq.write_text("")
    with pytest.raises(ValueError, match="No quality scores"):
        build_markov_model([str(fq)], max_reads=100, offset=33)