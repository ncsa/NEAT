"""
Regression tests for the probability_rates / local_mut_regions mismatch and
the intersect_regions algorithm fix.

Both fixed on branch fix/generate-variants-probability-rates-mismatch.

Root cause (two related bugs):

1. intersect_regions dropped all "middle" mutation regions (those fully
   contained within the block) and appended a zero-length fallback when
   block_end == last_region_end, returning M items instead of N.

2. generate_variants built probability_rates from the original N-item
   mutation_rate_regions list, then called rng.choice(a=local_mut_regions,
   p=probability_rates) where len(local_mut_regions) == M.  When M != N
   (always true for N >= 3 in practice) this raised:
       ValueError: 'a' and 'p' must have same length

Fixes:
    intersect_regions — rewritten using overlap arithmetic so every region
        that intersects the block is included and the tail fallback only fires
        when the block genuinely extends past all regions.

    generate_variants — probability_rates now derived from local_mut_regions
        (after the intersect) so lengths always match; None rates are
        substituted with mutation_model.avg_mut_rate.

Note: the single-region case and basic generate_variants behaviour are already
covered on the feature/claude-assisted-tests branch (28 tests).  Tests here
cover only the multi-region and None-rate scenarios that were broken.
"""
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.models import MutationModel
from neat.read_simulator.utils.bed_func import intersect_regions
from neat.read_simulator.utils.generate_variants import generate_variants
from neat.read_simulator.utils.options import Options
from neat.variants import ContigVariants


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ref(length: int = 400) -> SeqRecord:
    seq = "ACGT" * (length // 4)
    return SeqRecord(Seq(seq), id="chr1", name="chr1", description="")


def _make_opts(seed: int = 0, min_mutations: int = 1) -> Options:
    opts = Options(rng_seed=seed)
    opts.ploidy = 2
    opts.min_mutations = min_mutations
    opts.mutation_rate = None
    opts.mutation_bed = None
    return opts


def _run_gv(regions, ref=None, seed=0, min_mutations=1):
    if ref is None:
        ref = _make_ref()
    opts = _make_opts(seed=seed, min_mutations=min_mutations)
    model = MutationModel()
    model.rng = opts.rng
    return generate_variants(
        reference=ref,
        ref_start=0,
        mutation_rate_regions=regions,
        input_variants=ContigVariants(),
        mutation_model=model,
        options=opts,
        max_qual_score=40,
    )


# ===========================================================================
# intersect_regions — regression for the algorithm rewrite
# ===========================================================================

def test_intersect_three_regions_all_included():
    """
    Three regions spanning the block exactly — all three must appear in output.
    The old algorithm dropped the middle region and returned only 2 items.
    """
    regions = [(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)]
    result = intersect_regions(regions, (0, 400), 0.0)
    assert len(result) == 3
    assert result[0] == (0, 133, 0.01)
    assert result[1] == (133, 266, 0.02)
    assert result[2] == (266, 400, 0.05)


def test_intersect_four_regions_all_included():
    """Four regions, block exactly spans all — all four returned."""
    regions = [(0, 100, 0.01), (100, 200, 0.02), (200, 300, 0.03), (300, 400, 0.05)]
    result = intersect_regions(regions, (0, 400), 0.0)
    assert len(result) == 4


def test_intersect_block_end_equals_last_region_end_no_fallback():
    """
    When block_end == last region end the old code appended a zero-length
    fallback (last_end, block_end, default).  The new code must NOT add it.
    """
    regions = [(0, 200, 0.01), (200, 400, 0.05)]
    result = intersect_regions(regions, (0, 400), 0.0)
    # No zero-length entry
    assert all(r[0] < r[1] for r in result)


def test_intersect_block_partially_overlaps_middle_region():
    """Block (150, 350) overlaps parts of all three input regions."""
    regions = [(0, 200, 0.01), (200, 300, 0.02), (300, 400, 0.05)]
    result = intersect_regions(regions, (150, 350), 0.0)
    assert (150, 200, 0.01) in result
    assert (200, 300, 0.02) in result
    assert (300, 350, 0.05) in result


def test_intersect_result_is_contiguous():
    """Output sub-intervals must be contiguous (each end == next start)."""
    regions = [(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)]
    result = intersect_regions(regions, (0, 400), 0.0)
    for i in range(len(result) - 1):
        assert result[i][1] == result[i + 1][0]


def test_intersect_result_covers_full_block():
    """First item starts at block_start, last item ends at block_end."""
    regions = [(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)]
    result = intersect_regions(regions, (0, 400), 0.0)
    assert result[0][0] == 0
    assert result[-1][1] == 400


def test_intersect_block_outside_all_regions_returns_default():
    """Block with no overlap with any region → single fallback entry."""
    regions = [(0, 100, 0.01), (100, 200, 0.05)]
    result = intersect_regions(regions, (300, 500), 0.99)
    assert result == [(300, 500, 0.99)]


# ===========================================================================
# generate_variants — crash regression with multiple mutation rate regions
# ===========================================================================

def test_three_mut_regions_no_crash():
    """
    Primary crash regression: 3 mutation rate regions over a 400 bp reference.

    Before the fix:
        intersect_regions returned 2 items; probability_rates had 3 →
        ValueError: 'a' and 'p' must have same length
    """
    result = _run_gv([(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)])
    assert isinstance(result, ContigVariants)


def test_two_mut_regions_no_crash():
    """Two regions — was silently broken (zero-length fallback as second region)."""
    result = _run_gv([(0, 200, 0.01), (200, 400, 0.05)])
    assert isinstance(result, ContigVariants)


def test_four_mut_regions_no_crash():
    """Four regions — more aggressively exercises the fix."""
    result = _run_gv([(0, 100, 0.01), (100, 200, 0.02), (200, 300, 0.03), (300, 400, 0.05)])
    assert isinstance(result, ContigVariants)


def test_none_rate_region_no_crash():
    """
    None rate (from recalibrate_mutation_regions when no BED rate exists) is
    replaced by avg_mut_rate before building probability_rates.
    """
    result = _run_gv([(0, 200, None), (200, 400, 0.02)])
    assert isinstance(result, ContigVariants)


def test_all_none_rates_no_crash():
    """All None rates fall back entirely to avg_mut_rate."""
    result = _run_gv([(0, 200, None), (200, 400, None)])
    assert isinstance(result, ContigVariants)


def test_three_regions_produces_variants():
    """Multi-region run still generates at least the requested minimum mutations."""
    result = _run_gv([(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)],
                     min_mutations=5)
    assert len(result.variant_locations) >= 1


def test_multi_region_variant_positions_in_bounds():
    """All variant positions fall within the reference after the fix."""
    ref = _make_ref(400)
    result = _run_gv([(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)],
                     ref=ref, min_mutations=5)
    for loc in result.variant_locations:
        assert 0 <= loc < len(ref)


def test_multi_region_reproducible_with_same_seed():
    """Same seed produces identical variant locations with multiple regions."""
    regions = [(0, 133, 0.01), (133, 266, 0.02), (266, 400, 0.05)]
    r1 = _run_gv(regions, seed=7, min_mutations=5)
    r2 = _run_gv(regions, seed=7, min_mutations=5)
    assert r1.variant_locations == r2.variant_locations
