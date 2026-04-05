"""
Regression test for fix/contig-variants-remove-variant.

remove_variant was using variant.position instead of variant.position1,
causing it to silently no-op on any variant (since base variants only
define position1, not position).
"""
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.variants.contig_variants import ContigVariants
from neat.variants import SingleNucleotideVariant


_GT = np.array([0, 1])


def _snv(pos, alt="T"):
    return SingleNucleotideVariant(pos, alt, _GT.copy(), "42")


def test_remove_variant_removes_existing():
    cv = ContigVariants()
    v = _snv(10, "T")
    cv.add_variant(v)
    assert 10 in cv.variant_locations
    cv.remove_variant(v)
    assert 10 not in cv.variant_locations


def test_remove_variant_empties_contig_variants_list():
    cv = ContigVariants()
    v = _snv(5, "C")
    cv.add_variant(v)
    cv.remove_variant(v)
    # The location is removed from variant_locations; the dict entry is empty
    assert v not in cv.contig_variants.get(5, [])


def test_remove_variant_leaves_other_locations_intact():
    cv = ContigVariants()
    v1 = _snv(5, "C")
    v2 = _snv(20, "G")
    cv.add_variant(v1)
    cv.add_variant(v2)
    cv.remove_variant(v1)
    assert 20 in cv.variant_locations
    assert 5 not in cv.variant_locations


def test_remove_variant_with_multiple_at_same_location():
    """Removing one of two variants at the same position leaves the other."""
    cv = ContigVariants()
    v1 = SingleNucleotideVariant(10, "T", np.array([1, 0]), "42")
    v2 = SingleNucleotideVariant(10, "G", np.array([0, 1]), "42")
    cv.add_variant(v1)
    cv.add_variant(v2)
    cv.remove_variant(v1)
    assert 10 in cv.variant_locations
    assert v2 in cv.contig_variants[10]


def test_remove_variant_noop_when_absent():
    """Calling remove_variant on a variant not in ContigVariants is safe."""
    cv = ContigVariants()
    v = _snv(99, "A")
    cv.remove_variant(v)  # should not raise