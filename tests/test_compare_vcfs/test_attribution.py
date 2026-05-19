"""
Tests for neat/compare_vcfs/attribution.py — NEAT-aware FN classification.
"""
from pathlib import Path
from types import SimpleNamespace

import pytest

from neat.compare_vcfs.attribution import (
    REASON_OUTSIDE_CONTIGS,
    REASON_OUTSIDE_MUTATION_BED,
    REASON_OUTSIDE_TARGET_BED,
    REASON_UNKNOWN,
    attribute_fn,
    attribute_fns,
    load_bed_intervals,
    position_in_intervals,
)


# ---------------------------------------------------------------------------
# load_bed_intervals
# ---------------------------------------------------------------------------

def test_load_bed_intervals_returns_none_for_none_path():
    assert load_bed_intervals(None) is None


def test_load_bed_intervals_parses_basic_bed(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t100\t200\nchr1\t300\t400\nchr2\t50\t75\n")
    result = load_bed_intervals(bed)
    assert result == {"chr1": [(100, 200), (300, 400)], "chr2": [(50, 75)]}


def test_load_bed_intervals_skips_comments_and_headers(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text(
        "# header comment\n"
        "track name=foo\n"
        "browser position chr1\n"
        "\n"
        "chr1\t100\t200\n"
    )
    assert load_bed_intervals(bed) == {"chr1": [(100, 200)]}


def test_load_bed_intervals_skips_malformed_rows(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text(
        "chr1\t100\t200\n"
        "short_row\n"
        "chr1\tNaN\t300\n"
        "chr2\t10\t20\n"
    )
    assert load_bed_intervals(bed) == {"chr1": [(100, 200)], "chr2": [(10, 20)]}


def test_load_bed_intervals_sorts_intervals_per_contig(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t300\t400\nchr1\t100\t200\nchr1\t500\t600\n")
    result = load_bed_intervals(bed)
    assert result["chr1"] == [(100, 200), (300, 400), (500, 600)]


def test_load_bed_intervals_ignores_extra_columns(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t100\t200\tname1\t1000\t+\n")
    assert load_bed_intervals(bed) == {"chr1": [(100, 200)]}


def test_load_bed_intervals_returns_empty_dict_for_empty_file(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    assert load_bed_intervals(bed) == {}


# ---------------------------------------------------------------------------
# position_in_intervals — coordinate boundaries are the tricky part
# ---------------------------------------------------------------------------

def test_position_in_intervals_inside_interval():
    """VCF pos 150 → 0-based 149; falls in [100, 200)."""
    assert position_in_intervals(150, [(100, 200)]) is True


def test_position_in_intervals_at_start_inclusive():
    """BED start is inclusive. VCF pos 101 → 0-based 100; [100, 200) starts at 100."""
    assert position_in_intervals(101, [(100, 200)]) is True


def test_position_in_intervals_at_end_inclusive_for_1based():
    """BED end is exclusive. VCF pos 200 → 0-based 199; in [100, 200)."""
    assert position_in_intervals(200, [(100, 200)]) is True


def test_position_in_intervals_just_past_end():
    """VCF pos 201 → 0-based 200; [100, 200) does not include 200."""
    assert position_in_intervals(201, [(100, 200)]) is False


def test_position_in_intervals_before_first():
    assert position_in_intervals(50, [(100, 200)]) is False


def test_position_in_intervals_between_intervals():
    assert position_in_intervals(250, [(100, 200), (300, 400)]) is False


def test_position_in_intervals_picks_right_interval_in_sorted_list():
    intervals = [(100, 200), (300, 400), (500, 600)]
    assert position_in_intervals(350, intervals) is True
    assert position_in_intervals(550, intervals) is True
    assert position_in_intervals(450, intervals) is False


def test_position_in_intervals_handles_empty_list():
    assert position_in_intervals(100, []) is False


def test_position_in_intervals_finds_match_in_overlapping_intervals():
    """Overlapping intervals: ensure we don't miss containment in an earlier interval."""
    intervals = [(100, 500), (200, 300)]  # overlapping; sorted by start
    assert position_in_intervals(450, intervals) is True


# ---------------------------------------------------------------------------
# attribute_fn — single-record classification
# ---------------------------------------------------------------------------

def test_attribute_fn_returns_outside_contigs_when_chrom_not_simulated():
    reasons = attribute_fn("chrZ", 100, frozenset({"chr1"}), None, None)
    assert reasons == [REASON_OUTSIDE_CONTIGS]


def test_attribute_fn_outside_contigs_does_not_combine_with_bed_reasons():
    """If the contig isn't simulated, BED checks are skipped entirely."""
    reasons = attribute_fn("chrZ", 100, frozenset({"chr1"}),
                            {"chr1": [(0, 1000)]}, {"chr1": [(0, 1000)]})
    assert reasons == [REASON_OUTSIDE_CONTIGS]


def test_attribute_fn_unknown_when_no_beds_configured():
    reasons = attribute_fn("chr1", 100, frozenset({"chr1"}), None, None)
    assert reasons == [REASON_UNKNOWN]


def test_attribute_fn_outside_mutation_bed_only():
    reasons = attribute_fn(
        "chr1", 100, frozenset({"chr1"}),
        mutation_intervals={"chr1": [(500, 600)]},  # 100 outside
        target_intervals=None,
    )
    assert reasons == [REASON_OUTSIDE_MUTATION_BED]


def test_attribute_fn_outside_target_bed_only():
    reasons = attribute_fn(
        "chr1", 100, frozenset({"chr1"}),
        mutation_intervals=None,
        target_intervals={"chr1": [(500, 600)]},
    )
    assert reasons == [REASON_OUTSIDE_TARGET_BED]


def test_attribute_fn_multiple_reasons_combined():
    """A FN outside both beds gets both tags, in canonical order."""
    reasons = attribute_fn(
        "chr1", 100, frozenset({"chr1"}),
        mutation_intervals={"chr1": [(500, 600)]},
        target_intervals={"chr1": [(700, 800)]},
    )
    assert reasons == [REASON_OUTSIDE_MUTATION_BED, REASON_OUTSIDE_TARGET_BED]


def test_attribute_fn_unknown_when_inside_all_configured_beds():
    """If the FN is inside every configured bed, NEAT has no explanation."""
    reasons = attribute_fn(
        "chr1", 150, frozenset({"chr1"}),
        mutation_intervals={"chr1": [(100, 200)]},
        target_intervals={"chr1": [(100, 200)]},
    )
    assert reasons == [REASON_UNKNOWN]


def test_attribute_fn_bed_missing_chrom_is_outside():
    """A bed that doesn't mention this chrom counts as 'outside' for it."""
    reasons = attribute_fn(
        "chr1", 100, frozenset({"chr1"}),
        mutation_intervals={"chr2": [(0, 1000)]},  # chr1 absent
        target_intervals=None,
    )
    assert reasons == [REASON_OUTSIDE_MUTATION_BED]


# ---------------------------------------------------------------------------
# attribute_fns — the integration entry point used by the runner
# ---------------------------------------------------------------------------

def _fake_record(chrom: str, pos: int):
    return SimpleNamespace(chrom=chrom, pos=pos)


def test_attribute_fns_returns_one_tag_set_per_record(tmp_path):
    mut = tmp_path / "mut.bed"
    mut.write_text("chr1\t0\t1000\n")
    summary = {
        "delivered": {"contigs_simulated": ["chr1", "chr2"]},
        "config": {"mutation_bed": str(mut), "target_bed": None},
    }
    fns = [_fake_record("chr1", 500), _fake_record("chr2", 100), _fake_record("chrZ", 100)]
    result = attribute_fns(fns, summary)
    assert [r for _, r in result] == [
        [REASON_UNKNOWN],                  # chr1:500 inside mut bed
        [REASON_OUTSIDE_MUTATION_BED],     # chr2:100 — chr2 not in mut bed
        [REASON_OUTSIDE_CONTIGS],          # chrZ wasn't simulated
    ]


def test_attribute_fns_no_beds_configured(tmp_path):
    summary = {
        "delivered": {"contigs_simulated": ["chr1"]},
        "config": {"mutation_bed": None, "target_bed": None},
    }
    fns = [_fake_record("chr1", 100)]
    result = attribute_fns(fns, summary)
    assert result[0][1] == [REASON_UNKNOWN]


def test_attribute_fns_returns_empty_list_for_empty_input():
    summary = {"delivered": {"contigs_simulated": []}, "config": {}}
    assert attribute_fns([], summary) == []


def test_attribute_fns_pairs_each_record_with_reasons(tmp_path):
    summary = {
        "delivered": {"contigs_simulated": ["chr1"]},
        "config": {"mutation_bed": None, "target_bed": None},
    }
    rec = _fake_record("chr1", 42)
    [(returned_rec, reasons)] = attribute_fns([rec], summary)
    assert returned_rec is rec
    assert reasons == [REASON_UNKNOWN]
