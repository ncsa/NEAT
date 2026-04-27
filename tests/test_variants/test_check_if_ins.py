"""
Regression test for fix/contig-variants-check-if-ins.

check_if_ins was passing the variant object to Insertion.contains() which
expects an int (a position). This caused the method to always return None
even when a variant's position falls within an insertion's span.
"""
import numpy as np

from neat.variants.contig_variants import ContigVariants
from neat.variants import Insertion, SingleNucleotideVariant


_GT = np.array([0, 1])


def _ins(pos, alt="ACGTT", length=4, gt=None):
    return Insertion(pos, length, alt, gt if gt is not None else _GT.copy(), "42")


def _snv(pos, alt="T", gt=None):
    return SingleNucleotideVariant(pos, alt, gt if gt is not None else _GT.copy(), "42")


def test_check_if_ins_finds_containing_insertion():
    """SNV inside insertion span is detected."""
    cv = ContigVariants()
    ins = _ins(10, "ACGTT", 4)
    cv.add_variant(ins)
    snv = _snv(11)  # position1=11 is within [10, 14)
    result = cv.check_if_ins(snv)
    assert result is ins


def test_check_if_ins_at_insertion_start():
    """SNV at the insertion's own position is detected."""
    cv = ContigVariants()
    ins = _ins(10, "ACGTT", 4)
    cv.add_variant(ins)
    snv = _snv(10)
    assert cv.check_if_ins(snv) is ins


def test_check_if_ins_at_insertion_end_exclusive():
    """Position at insertion start + length is outside the span."""
    cv = ContigVariants()
    ins = _ins(10, "ACGTT", 4)  # spans [10, 14)
    cv.add_variant(ins)
    snv = _snv(14)
    assert cv.check_if_ins(snv) is None


def test_check_if_ins_outside_span_returns_none():
    cv = ContigVariants()
    ins = _ins(10, "ACGTT", 4)
    cv.add_variant(ins)
    snv = _snv(50)
    assert cv.check_if_ins(snv) is None


def test_check_if_ins_empty_returns_none():
    cv = ContigVariants()
    assert cv.check_if_ins(_snv(10)) is None


def test_check_if_ins_genotype_mismatch_returns_none():
    """Even if position matches, different genotype means no match."""
    cv = ContigVariants()
    ins = _ins(10, "ACGTT", 4, gt=np.array([1, 0]))
    cv.add_variant(ins)
    snv = _snv(11, gt=np.array([0, 1]))  # different genotype
    assert cv.check_if_ins(snv) is None
