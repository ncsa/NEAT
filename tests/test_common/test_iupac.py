"""
Tests for IUPAC ambiguity-code resolution (issue #291).
"""
import numpy as np
from numpy.random import default_rng

from neat.common.constants_and_defaults import IUPAC_CODES, resolve_iupac_bases


def test_no_ambiguity_codes_returns_unchanged():
    """A plain A/C/G/T/N sequence is returned untouched with a zero count."""
    seq = "ACGTACGTNNNNACGT"
    out, n = resolve_iupac_bases(seq, default_rng(0))
    assert out == seq
    assert n == 0


def test_n_is_not_resolved():
    """'N' is deliberately left alone (it has its own downstream masking)."""
    seq = "N" * 50
    out, n = resolve_iupac_bases(seq, default_rng(1))
    assert out == seq
    assert n == 0


def test_every_code_resolves_into_its_allowed_set():
    """Each ambiguity code is replaced by one of the bases it represents."""
    rng = default_rng(7)
    for code, bases in IUPAC_CODES.items():
        seq = code * 200
        out, n = resolve_iupac_bases(seq, rng)
        assert n == 200
        assert set(out) <= set(bases), f"{code} produced bases outside {bases}: {set(out)}"
        # Sanity: nothing left unresolved.
        assert code not in out


def test_count_matches_number_of_codes():
    """The returned count equals the number of ambiguity positions resolved."""
    seq = "ACGTRYSWKMACGTBDHVN"  # 10 ambiguity codes embedded
    out, n = resolve_iupac_bases(seq, default_rng(3))
    assert n == 10
    # Only A/C/G/T/N remain.
    assert set(out) <= set("ACGTN")
    # The N is preserved in place (last char).
    assert out[-1] == "N"


def test_resolution_is_reproducible_with_seed():
    """Same seed → identical resolution."""
    seq = "RYSWKMBDHV" * 100
    a, na = resolve_iupac_bases(seq, default_rng(99))
    b, nb = resolve_iupac_bases(seq, default_rng(99))
    assert a == b
    assert na == nb == 1000


def test_mixed_distribution_covers_both_alternatives():
    """Over many positions, a 2-base code resolves to both of its options."""
    seq = "R" * 1000
    out, _ = resolve_iupac_bases(seq, default_rng(2024))
    assert set(out) == {"A", "G"}, f"expected both A and G, got {set(out)}"


def test_length_preserved():
    """Resolution never changes sequence length."""
    seq = "ACGTRYSWKMNBDHV"
    out, _ = resolve_iupac_bases(seq, default_rng(5))
    assert len(out) == len(seq)


def test_lowercase_not_resolved():
    """Resolution operates on uppercased input; lowercase ambiguity codes pass through.

    The caller (split_inputs) uppercases before calling, so this just documents that the
    function itself only acts on the canonical uppercase codes and won't silently mangle
    unexpected input.
    """
    seq = "rystACGT"
    out, n = resolve_iupac_bases(seq, default_rng(0))
    assert out == seq
    assert n == 0
