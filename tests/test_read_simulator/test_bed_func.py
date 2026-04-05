"""
Tests for neat/read_simulator/utils/bed_func.py
"""
import pytest
from pathlib import Path

from neat.read_simulator.utils.bed_func import (
    intersect_regions,
    recalibrate_mutation_regions,
    fill_out_bed_dict,
    fill_out_mut_regions,
    parse_single_bed,
    parse_beds,
)
from neat.read_simulator.utils.options import Options


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_REF = {"chr1": 1000, "chr2": 500}

_TARGET_KEY  = ("target",  True,  True)
_DISCARD_KEY = ("discard", False, True)
_MUT_KEY     = ("mutation", None, True)
_NO_BED_TARGET  = ("target",  True,  False)
_NO_BED_DISCARD = ("discard", False, False)


def _write_bed(tmp_path: Path, name: str, lines: list[str]) -> Path:
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


# ===========================================================================
# intersect_regions
# ===========================================================================

def test_intersect_block_within_single_region():
    regions = [(0, 1000, 0.01)]
    result = intersect_regions(regions, (100, 200), 0.0)
    assert result == [(100, 200, 0.01)]


def test_intersect_block_exactly_fills_region():
    regions = [(0, 500, 0.02), (500, 1000, 0.05)]
    result = intersect_regions(regions, (0, 500), 0.0)
    assert result == [(0, 500, 0.02)]


def test_intersect_block_spans_two_regions():
    regions = [(0, 300, 0.01), (300, 700, 0.05), (700, 1000, 0.02)]
    result = intersect_regions(regions, (200, 400), 0.0)
    assert result == [(200, 300, 0.01), (300, 400, 0.05)]


def test_intersect_block_start_aligns_with_region_boundary():
    regions = [(0, 500, 0.01), (500, 1000, 0.03)]
    result = intersect_regions(regions, (500, 800), 0.0)
    assert result == [(500, 800, 0.03)]


def test_intersect_block_extends_beyond_all_regions():
    regions = [(0, 500, 0.01), (500, 800, 0.02)]
    result = intersect_regions(regions, (600, 1000), 0.0)
    # Ends past last region — appends fallback using default_value
    assert result[-1][2] == 0.0
    assert result[-1][1] == 1000


# ===========================================================================
# fill_out_mut_regions
# ===========================================================================

def test_fill_out_mut_regions_empty():
    result = fill_out_mut_regions([], (0, 1000), 0.01)
    assert result == [(0, 1000, 0.01)]


def test_fill_out_mut_regions_full_coverage():
    regions = [(0, 1000, 0.05)]
    result = fill_out_mut_regions(regions, (0, 1000), 0.01)
    assert result == [(0, 1000, 0.05)]


def test_fill_out_mut_regions_gap_before():
    regions = [(200, 800, 0.05)]
    result = fill_out_mut_regions(regions, (0, 1000), 0.01)
    assert result[0] == (0, 200, 0.01)
    assert result[1] == (200, 800, 0.05)
    assert result[2] == (800, 1000, 0.01)


def test_fill_out_mut_regions_gap_after():
    regions = [(0, 600, 0.05)]
    result = fill_out_mut_regions(regions, (0, 1000), 0.01)
    assert result[-1] == (600, 1000, 0.01)


def test_fill_out_mut_regions_multiple_gaps():
    regions = [(100, 300, 0.01), (600, 800, 0.02)]
    result = fill_out_mut_regions(regions, (0, 1000), 0.0)
    starts = [r[0] for r in result]
    ends   = [r[1] for r in result]
    # Regions must be contiguous from 0 to 1000
    assert starts[0] == 0
    assert ends[-1] == 1000
    for i in range(len(result) - 1):
        assert result[i][1] == result[i + 1][0]


# ===========================================================================
# recalibrate_mutation_regions
# ===========================================================================

def test_recalibrate_none_values_replaced():
    regions = [(0, 500, None), (500, 1000, 0.02)]
    result = recalibrate_mutation_regions(regions, (0, 1000), 0.01)
    rates = [r[2] for r in result]
    assert None not in rates


def test_recalibrate_block_contained_in_one_region():
    regions = [(0, 1000, 0.05)]
    result = recalibrate_mutation_regions(regions, (200, 400), 0.01)
    assert result == [(200, 400, 0.05)]


def test_recalibrate_block_spans_two_regions():
    regions = [(0, 500, 0.01), (500, 1000, 0.03)]
    result = recalibrate_mutation_regions(regions, (300, 700), 0.0)
    assert any(r[2] == 0.01 for r in result)
    assert any(r[2] == 0.03 for r in result)


def test_recalibrate_result_is_contiguous():
    regions = [(0, 400, 0.01), (400, 600, None), (600, 1000, 0.02)]
    result = recalibrate_mutation_regions(regions, (0, 1000), 0.005)
    for i in range(len(result) - 1):
        assert result[i][1] == result[i + 1][0], "Gap between regions"


def test_recalibrate_result_covers_full_block():
    regions = [(0, 1000, 0.01)]
    coords = (100, 900)
    result = recalibrate_mutation_regions(regions, coords, 0.005)
    assert result[0][0] == 100
    assert result[-1][1] == 900


# ===========================================================================
# parse_single_bed — no bed file (uniform coverage)
# ===========================================================================

def test_parse_single_bed_no_file_target():
    result = parse_single_bed(None, _REF, _NO_BED_TARGET)
    assert set(result.keys()) == {"chr1", "chr2"}
    assert result["chr1"] == [(0, 1000, True)]
    assert result["chr2"] == [(0, 500, True)]


def test_parse_single_bed_no_file_discard():
    result = parse_single_bed(None, _REF, _NO_BED_DISCARD)
    assert result["chr1"] == [(0, 1000, False)]


# ===========================================================================
# parse_single_bed — target BED
# ===========================================================================

def test_parse_single_bed_target(tmp_path):
    bed = _write_bed(tmp_path, "target.bed", [
        "chr1\t100\t400",
        "chr1\t600\t800",
    ])
    result = parse_single_bed(str(bed), _REF, _TARGET_KEY)
    entries = result["chr1"]
    assert len(entries) == 2
    assert entries[0] == (100, 400, True)
    assert entries[1] == (600, 800, True)


def test_parse_single_bed_target_skips_comments(tmp_path):
    bed = _write_bed(tmp_path, "target.bed", [
        "# comment line",
        "chr1\t0\t500",
    ])
    result = parse_single_bed(str(bed), _REF, _TARGET_KEY)
    assert result["chr1"] == [(0, 500, True)]


def test_parse_single_bed_target_chrom_not_in_ref(tmp_path):
    bed = _write_bed(tmp_path, "target.bed", [
        "chrX\t0\t500",
        "chr1\t0\t500",
    ])
    result = parse_single_bed(str(bed), _REF, _TARGET_KEY)
    # chrX is skipped, chr1 is kept
    assert result["chr1"] == [(0, 500, True)]
    assert "chrX" not in result


# ===========================================================================
# parse_single_bed — discard BED
# ===========================================================================

def test_parse_single_bed_discard(tmp_path):
    bed = _write_bed(tmp_path, "discard.bed", ["chr1\t200\t600"])
    result = parse_single_bed(str(bed), _REF, _DISCARD_KEY)
    assert result["chr1"] == [(200, 600, True)]


# ===========================================================================
# parse_single_bed — mutation BED
# ===========================================================================

def test_parse_single_bed_mutation(tmp_path):
    bed = _write_bed(tmp_path, "mut.bed", [
        "chr1\t0\t500\tmut_rate=0.001",
        "chr1\t500\t1000\tmut_rate=0.005",
    ])
    result = parse_single_bed(str(bed), _REF, _MUT_KEY)
    assert result["chr1"][0] == (0, 500, 0.001)
    assert result["chr1"][1] == (500, 1000, 0.005)


def test_parse_single_bed_mutation_extra_metadata(tmp_path):
    """mut_rate can appear among other semicolon-delimited metadata."""
    bed = _write_bed(tmp_path, "mut.bed", [
        "chr1\t0\t500\tfoo=bar;mut_rate=0.002;baz=qux",
    ])
    result = parse_single_bed(str(bed), _REF, _MUT_KEY)
    assert result["chr1"][0][2] == pytest.approx(0.002)


def test_parse_single_bed_mutation_high_rate_logs_warning(tmp_path, caplog):
    """Mutation rate > 0.3 triggers a warning log (bed_func.py line 209)."""
    import logging
    bed = _write_bed(tmp_path, "mut.bed", ["chr1\t0\t500\tmut_rate=0.4"])
    with caplog.at_level(logging.WARNING, logger="neat.read_simulator.utils.bed_func"):
        result = parse_single_bed(str(bed), _REF, _MUT_KEY)
    assert "0.3" in caplog.text or "unusual" in caplog.text.lower()
    assert result["chr1"][0][2] == pytest.approx(0.4)


def test_parse_single_bed_mutation_missing_mut_rate_exits(tmp_path):
    bed = _write_bed(tmp_path, "mut.bed", ["chr1\t0\t500\tno_rate_here"])
    with pytest.raises(SystemExit):
        parse_single_bed(str(bed), _REF, _MUT_KEY)


def test_parse_single_bed_mutation_invalid_rate_value_exits(tmp_path):
    bed = _write_bed(tmp_path, "mut.bed", ["chr1\t0\t500\tmut_rate=notanumber"])
    with pytest.raises(SystemExit):
        parse_single_bed(str(bed), _REF, _MUT_KEY)


def test_parse_single_bed_malformed_line_exits(tmp_path):
    bed = _write_bed(tmp_path, "bad.bed", ["only_one_column"])
    with pytest.raises(SystemExit):
        parse_single_bed(str(bed), _REF, _TARGET_KEY)


# ===========================================================================
# fill_out_bed_dict
# ===========================================================================

def test_fill_out_bed_dict_empty_regions():
    region_dict = {"chr1": [], "chr2": []}
    result = fill_out_bed_dict(_REF, region_dict, _TARGET_KEY)
    assert result["chr1"] == [(0, 1000, True)]
    assert result["chr2"] == [(0, 500, True)]


def test_fill_out_bed_dict_target_fills_gaps_with_false():
    region_dict = {"chr1": [(200, 600, True)], "chr2": []}
    result = fill_out_bed_dict(_REF, region_dict, _TARGET_KEY)
    regions = result["chr1"]
    # Gap before (0–200) should be False; target region True; gap after (600–1000) False
    assert regions[0] == (0, 200, False)
    assert regions[1] == (200, 600, True)
    assert regions[2] == (600, 1000, False)


def test_fill_out_bed_dict_discard_fills_gaps_with_default():
    region_dict = {"chr1": [(200, 600, True)], "chr2": []}
    result = fill_out_bed_dict(_REF, region_dict, _DISCARD_KEY)
    regions = result["chr1"]
    assert regions[0] == (0, 200, False)   # not discarded
    assert regions[1] == (200, 600, True)  # discarded
    assert regions[2] == (600, 1000, False)


def test_fill_out_bed_dict_region_starts_at_zero():
    region_dict = {"chr1": [(0, 500, True)], "chr2": []}
    result = fill_out_bed_dict(_REF, region_dict, _TARGET_KEY)
    regions = result["chr1"]
    assert regions[0] == (0, 500, True)
    assert regions[1] == (500, 1000, False)  # gap filled


def test_fill_out_bed_dict_region_ends_at_contig_end():
    region_dict = {"chr1": [(500, 1000, True)], "chr2": []}
    result = fill_out_bed_dict(_REF, region_dict, _TARGET_KEY)
    regions = result["chr1"]
    assert regions[0] == (0, 500, False)
    assert regions[1] == (500, 1000, True)


def test_fill_out_bed_dict_contiguous():
    """Every filled dict should have contiguous, non-overlapping regions."""
    region_dict = {"chr1": [(100, 400, True), (600, 900, True)], "chr2": []}
    result = fill_out_bed_dict(_REF, region_dict, _TARGET_KEY)
    regions = result["chr1"]
    assert regions[0][0] == 0
    assert regions[-1][1] == 1000
    for i in range(len(regions) - 1):
        assert regions[i][1] == regions[i + 1][0]


def test_fill_out_bed_dict_multiple_contigs():
    region_dict = {
        "chr1": [(0, 500, True)],
        "chr2": [(100, 300, True)],
    }
    result = fill_out_bed_dict(_REF, region_dict, _TARGET_KEY)
    assert result["chr1"][-1][1] == 1000
    assert result["chr2"][-1][1] == 500


# ===========================================================================
# parse_beds — high-level integration of the above
# ===========================================================================

def test_parse_beds_no_beds():
    opts = Options(rng_seed=0)
    opts.target_bed = None
    opts.discard_bed = None
    opts.mutation_bed = None
    target_dict, discard_dict, mut_dict = parse_beds(opts, _REF)
    # With no beds: target is all-True, discard is all-False, mut is all-None
    assert target_dict["chr1"] == [(0, 1000, True)]
    assert discard_dict["chr1"] == [(0, 1000, False)]
    assert mut_dict["chr1"] == [(0, 1000, None)]


def test_parse_beds_with_target_bed(tmp_path):
    bed = _write_bed(tmp_path, "target.bed", ["chr1\t200\t800"])
    opts = Options(rng_seed=0)
    opts.target_bed = str(bed)
    opts.discard_bed = None
    opts.mutation_bed = None
    target_dict, discard_dict, mut_dict = parse_beds(opts, _REF)
    # Target regions outside bed are filled with False
    target_chr1 = target_dict["chr1"]
    true_regions = [r for r in target_chr1 if r[2] is True]
    assert true_regions == [(200, 800, True)]


def test_parse_beds_with_discard_bed(tmp_path):
    bed = _write_bed(tmp_path, "discard.bed", ["chr1\t300\t700"])
    opts = Options(rng_seed=0)
    opts.target_bed = None
    opts.discard_bed = str(bed)
    opts.mutation_bed = None
    target_dict, discard_dict, mut_dict = parse_beds(opts, _REF)
    discard_chr1 = discard_dict["chr1"]
    discarded = [r for r in discard_chr1 if r[2] is True]
    assert discarded == [(300, 700, True)]


def test_parse_beds_with_mutation_bed(tmp_path):
    bed = _write_bed(tmp_path, "mut.bed", [
        "chr1\t0\t500\tmut_rate=0.01",
        "chr1\t500\t1000\tmut_rate=0.02",
    ])
    opts = Options(rng_seed=0)
    opts.target_bed = None
    opts.discard_bed = None
    opts.mutation_bed = str(bed)
    target_dict, discard_dict, mut_dict = parse_beds(opts, _REF)
    assert mut_dict["chr1"][0][2] == pytest.approx(0.01)
    assert mut_dict["chr1"][1][2] == pytest.approx(0.02)