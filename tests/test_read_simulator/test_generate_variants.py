"""
Unit tests for neat/read_simulator/utils/generate_variants.py
"""
import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.models import MutationModel
from neat.read_simulator.utils.generate_variants import (
    find_random_non_n,
    map_non_n_regions,
    generate_variants,
)
from neat.read_simulator.utils.options import Options
from neat.variants import ContigVariants, SingleNucleotideVariant, Insertion, Deletion


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CLEAN_SEQ = "ACGT" * 50          # 200 bp, no N's
_N_HEAVY   = "N" * 95 + "ACGT"    # 99 bases, >10% N

def _make_reference(seq: str = _CLEAN_SEQ, name: str = "chr1") -> SeqRecord:
    return SeqRecord(Seq(seq), id=name, name=name, description="")


def _make_options(seed: int = 42, mutation_rate: float = None) -> Options:
    opts = Options(rng_seed=seed)
    opts.ploidy = 2
    opts.mutation_rate = mutation_rate
    opts.min_mutations = 1
    return opts


def _make_model() -> MutationModel:
    return MutationModel()


def _full_rate_regions(seq_len: int, rate: float = 0.01, offset: int = 0):
    """Single mutation-rate region covering [offset, offset+seq_len)."""
    return [(offset, offset + seq_len, rate)]


# ===========================================================================
# find_random_non_n
# ===========================================================================

def test_find_random_non_n_returns_valid_index():
    rng = np.random.default_rng(0)
    safe_zones = np.ones(10, dtype=int)
    idx = find_random_non_n(rng, safe_zones)
    assert 0 <= idx < 10


def test_find_random_non_n_respects_probabilities():
    """With all weight on position 3, should always return 3."""
    rng = np.random.default_rng(0)
    safe_zones = np.zeros(10, dtype=int)
    safe_zones[3] = 1
    for _ in range(20):
        assert find_random_non_n(rng, safe_zones) == 3


def test_find_random_non_n_reproducible_with_seed():
    safe_zones = np.array([1, 2, 3, 4, 5], dtype=float)
    r1 = find_random_non_n(np.random.default_rng(7), safe_zones.copy())
    r2 = find_random_non_n(np.random.default_rng(7), safe_zones.copy())
    assert r1 == r2


def test_find_random_non_n_single_element():
    rng = np.random.default_rng(0)
    safe_zones = np.array([5], dtype=float)
    assert find_random_non_n(rng, safe_zones) == 0


# ===========================================================================
# map_non_n_regions
# ===========================================================================

def test_map_non_n_regions_clean_sequence():
    result = map_non_n_regions("ACGTACGTACGT")
    assert len(result) == 12
    assert all(result == 1)


def test_map_non_n_regions_single_n():
    # 1 N in 50 bases = 2% N → valid map returned
    seq = "A" * 24 + "N" + "A" * 25
    result = map_non_n_regions(seq)
    assert len(result) == 50
    assert result[24] == 0
    assert result[0] == 1
    assert result[49] == 1


def test_map_non_n_regions_too_many_ns_returns_empty():
    """More than 10% N → returns empty array."""
    seq = "N" * 50 + "ACGT" * 5  # 70 bp, 71% N
    result = map_non_n_regions(seq)
    assert len(result) == 0


def test_map_non_n_regions_exactly_at_threshold():
    """Exactly 10% N → should return empty (condition is <= 0.90 non-N)."""
    # 10 N's + 90 ACGT → 90% non-N, average == 0.90, which hits the <= boundary
    seq = "N" * 10 + "A" * 90
    result = map_non_n_regions(seq)
    assert len(result) == 0


def test_map_non_n_regions_all_n_returns_empty():
    result = map_non_n_regions("NNNNNNNNNN")
    assert len(result) == 0


def test_map_non_n_regions_run_of_ns():
    seq = "ACGT" + "NNNN" + "ACGT"  # 12 bp, 4/12 ≈ 33% N → empty
    result = map_non_n_regions(seq)
    assert len(result) == 0


def test_map_non_n_regions_short_n_run_in_long_sequence():
    """2 N's in 100 bp (2% N) → valid map returned."""
    seq = "A" * 49 + "NN" + "A" * 49
    result = map_non_n_regions(seq)
    assert len(result) == 100
    assert result[49] == 0
    assert result[50] == 0
    assert result[0] == 1


# ===========================================================================
# generate_variants — input variants are copied into output
# ===========================================================================

def test_generate_variants_returns_contig_variants():
    ref = _make_reference()
    model = _make_model()
    opts = _make_options()
    opts.min_mutations = 0
    iv = ContigVariants()
    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), iv, model, opts, 40)
    assert isinstance(result, ContigVariants)


def test_generate_variants_input_variant_in_range_is_included():
    ref = _make_reference()
    model = _make_model()
    opts = _make_options()
    opts.min_mutations = 0

    iv = ContigVariants()
    snv = SingleNucleotideVariant(10, "T", np.array([0, 1]), "42")
    iv.add_variant(snv)

    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), iv, model, opts, 40)
    assert 10 in result


def test_generate_variants_input_variant_outside_range_excluded():
    ref = _make_reference(_CLEAN_SEQ[:100])
    model = _make_model()
    opts = _make_options()
    opts.min_mutations = 0

    iv = ContigVariants()
    # Place variant at position 150, but ref is only 100 bp starting at 0
    snv = SingleNucleotideVariant(150, "T", np.array([0, 1]), "42")
    iv.add_variant(snv)

    result = generate_variants(ref, 0, _full_rate_regions(100), iv, model, opts, 40)
    assert 150 not in result


def test_generate_variants_input_variant_offset_range():
    """ref_start=100, variant at 110 → should be included (100 <= 110 < 300)."""
    ref = _make_reference()  # 200 bp
    model = _make_model()
    opts = _make_options()
    opts.min_mutations = 0

    iv = ContigVariants()
    snv = SingleNucleotideVariant(110, "T", np.array([0, 1]), "42")
    iv.add_variant(snv)

    result = generate_variants(ref, 100, _full_rate_regions(200, offset=100), iv, model, opts, 40)
    assert 110 in result


def test_generate_variants_input_variant_before_offset_excluded():
    """ref_start=100, variant at 50 → should not be included."""
    ref = _make_reference()  # 200 bp
    model = _make_model()
    opts = _make_options()
    opts.min_mutations = 0

    iv = ContigVariants()
    snv = SingleNucleotideVariant(50, "T", np.array([0, 1]), "42")
    iv.add_variant(snv)

    result = generate_variants(ref, 100, _full_rate_regions(200, offset=100), iv, model, opts, 40)
    assert 50 not in result


# ===========================================================================
# generate_variants — random mutation generation
# ===========================================================================

def test_generate_variants_adds_at_least_min_mutations():
    """With min_mutations=1, at least 1 variant should be added."""
    ref = _make_reference()
    model = _make_model()
    opts = _make_options(seed=0)
    opts.min_mutations = 1

    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), ContigVariants(), model, opts, 40)
    assert len(result.variant_locations) >= 1


def test_generate_variants_qual_score_applied():
    """Randomly generated variants should have qual score == max_qual_score."""
    ref = _make_reference()
    model = _make_model()
    opts = _make_options(seed=1)
    opts.min_mutations = 5

    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), ContigVariants(), model, opts, 60)
    for loc in result.variant_locations:
        for var in result.contig_variants[loc]:
            # Input variants may keep their original score; randomly generated ones get max_qual_score
            if not getattr(var, 'is_input', False):
                assert var.qual_score == 60


def test_generate_variants_zero_mutation_rate_still_runs():
    """A mutation rate of 0 should not crash."""
    ref = _make_reference()
    model = _make_model()
    opts = _make_options(seed=2)
    opts.min_mutations = 0

    result = generate_variants(ref, 0, [(0, len(ref), 0.0)], ContigVariants(), model, opts, 40)
    assert isinstance(result, ContigVariants)


def test_generate_variants_high_mutation_rate_adds_many():
    """A high mutation rate should produce more variants than a low one."""
    ref = _make_reference("ACGT" * 100)  # 400 bp
    model_low = MutationModel(avg_mut_rate=0.001)
    model_high = MutationModel(avg_mut_rate=0.05)
    opts_low  = _make_options(seed=5)
    opts_high = _make_options(seed=5)
    opts_low.min_mutations  = 0
    opts_high.min_mutations = 0

    result_low  = generate_variants(ref, 0, _full_rate_regions(400, 0.001), ContigVariants(), model_low,  opts_low,  40)
    result_high = generate_variants(ref, 0, _full_rate_regions(400, 0.05),  ContigVariants(), model_high, opts_high, 40)
    assert len(result_high.variant_locations) >= len(result_low.variant_locations)


def test_generate_variants_variants_within_reference_bounds():
    """All generated variant positions should fall within [ref_start, ref_start + len(ref))."""
    seq = "ACGT" * 75  # 300 bp
    ref = _make_reference(seq)
    model = _make_model()
    opts = _make_options(seed=3)
    opts.min_mutations = 10

    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), ContigVariants(), model, opts, 40)
    for loc in result.variant_locations:
        assert 0 <= loc < len(seq), f"Variant at {loc} is out of range [0, {len(seq)})"


def test_generate_variants_reproducible_with_same_seed():
    ref = _make_reference()
    model1 = _make_model()
    model2 = _make_model()

    opts1 = _make_options(seed=99)
    opts2 = _make_options(seed=99)
    opts1.min_mutations = opts2.min_mutations = 5

    r1 = generate_variants(ref, 0, _full_rate_regions(len(ref)), ContigVariants(), model1, opts1, 40)
    r2 = generate_variants(ref, 0, _full_rate_regions(len(ref)), ContigVariants(), model2, opts2, 40)
    assert r1.variant_locations == r2.variant_locations


def test_generate_variants_variant_types_are_valid():
    """Every generated variant should be Insertion, Deletion, or SNV."""
    ref = _make_reference()
    model = _make_model()
    opts = _make_options(seed=7)
    opts.min_mutations = 10

    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), ContigVariants(), model, opts, 40)
    for loc in result.variant_locations:
        for var in result.contig_variants[loc]:
            assert isinstance(var, (Insertion, Deletion, SingleNucleotideVariant))


def test_generate_variants_single_rate_region_at_offset():
    """generate_variants works correctly when ref_start places the block mid-reference."""
    seq = "ACGT" * 50   # 200 bp
    ref = _make_reference(seq)
    model = _make_model()
    opts = _make_options(seed=11)
    opts.min_mutations = 3
    # Single rate region that exactly spans the offset block [200, 400)
    rate_regions = [(200, 400, 0.01)]

    result = generate_variants(ref, 200, rate_regions, ContigVariants(), model, opts, 40)
    assert isinstance(result, ContigVariants)
    assert len(result.variant_locations) >= 1


def test_generate_variants_input_and_random_together():
    """Pre-existing input variant + randomly generated variants all appear in output."""
    ref = _make_reference()
    model = _make_model()
    opts = _make_options(seed=13)
    opts.min_mutations = 3

    iv = ContigVariants()
    snv = SingleNucleotideVariant(5, "C", np.array([1, 0]), "42")
    iv.add_variant(snv)

    result = generate_variants(ref, 0, _full_rate_regions(len(ref)), iv, model, opts, 40)
    # Input variant must still be present
    assert 5 in result
    # And at least some random variants were added
    assert len(result.variant_locations) > 1
