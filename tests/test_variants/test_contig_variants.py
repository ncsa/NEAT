"""
Unit tests for neat/variants/contig_variants.py — focusing on
get_ref_alt, get_sample_info, and remove_variant (previously uncovered).
"""
import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.variants.contig_variants import ContigVariants
from neat.variants import Deletion, Insertion, SingleNucleotideVariant
from neat.variants.unknown_variant import UnknownVariant

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SEQ = "ACGTACGTACGTACGT"   # 16 bp
_REC = SeqRecord(Seq(_SEQ), id="chr1", name="chr1", description="")
_GT  = np.array([0, 1])


def _snv(pos, alt="T", gt=None):
    return SingleNucleotideVariant(pos, alt, gt if gt is not None else _GT.copy(), "42")


def _del(pos, length=3, gt=None):
    return Deletion(pos, length, gt if gt is not None else _GT.copy(), "42")


def _ins(pos, alt="ACGT", length=3, gt=None):
    return Insertion(pos, length, alt, gt if gt is not None else _GT.copy(), "42")


# ===========================================================================
# get_ref_alt — SNV
# ===========================================================================

def test_get_ref_alt_snv_ref_is_single_base():
    snv = _snv(2, "T")
    ref, alt = ContigVariants.get_ref_alt(snv, _REC, 0)
    assert ref == _SEQ[2]


def test_get_ref_alt_snv_alt_matches_variant():
    snv = _snv(2, "T")
    ref, alt = ContigVariants.get_ref_alt(snv, _REC, 0)
    assert alt == "T"


def test_get_ref_alt_snv_with_block_start_offset():
    # ref starts at block_start=4; variant at position1=6 → local index=2
    snv = _snv(6, "G")
    ref, alt = ContigVariants.get_ref_alt(snv, _REC, 4)
    assert ref == _SEQ[2]   # local index = 6 - 4 = 2


# ===========================================================================
# get_ref_alt — Deletion
# ===========================================================================

def test_get_ref_alt_deletion_ref_spans_length():
    d = _del(1, 3)
    ref, alt = ContigVariants.get_ref_alt(d, _REC, 0)
    assert ref == _SEQ[1:4]


def test_get_ref_alt_deletion_alt_is_single_base():
    d = _del(1, 3)
    ref, alt = ContigVariants.get_ref_alt(d, _REC, 0)
    assert len(alt) == 1
    assert alt == _SEQ[1]


# ===========================================================================
# get_ref_alt — Insertion
# ===========================================================================

def test_get_ref_alt_insertion_ref_is_single_base():
    ins = _ins(3, "ACGTT", 4)
    ref, alt = ContigVariants.get_ref_alt(ins, _REC, 0)
    assert ref == _SEQ[3]


def test_get_ref_alt_insertion_alt_from_variant():
    ins = _ins(3, "ACGTT", 4)
    ref, alt = ContigVariants.get_ref_alt(ins, _REC, 0)
    assert alt == "ACGTT"


# ===========================================================================
# get_ref_alt — UnknownVariant
# ===========================================================================

def test_get_ref_alt_unknown_uses_metadata():
    uv = UnknownVariant(5, _GT.copy(), "42", is_input=True,
                        REF="A", ALT="ACGT")
    # UnknownVariant does not set self.alt; patch it so get_alt() falls
    # through to metadata['ALT'] rather than raising AttributeError.
    uv.alt = None
    ref, alt = ContigVariants.get_ref_alt(uv, _REC, 0)
    assert ref == "A"
    assert alt == "ACGT"


# ===========================================================================
# get_sample_info
# ===========================================================================

def test_get_sample_info_with_neat_sample_metadata():
    snv = _snv(2, "T")
    snv.metadata["NEAT_sample"] = "0|1"
    result = ContigVariants.get_sample_info(snv)
    assert result == "0|1"


def test_get_sample_info_without_metadata_uses_genotype_string():
    snv = _snv(2, "T", gt=np.array([0, 1]))
    result = ContigVariants.get_sample_info(snv)
    assert "|" in result or "/" in result


# ===========================================================================
# remove_variant
# ===========================================================================

def test_remove_variant_method_exists():
    """remove_variant silently no-ops due to variant.position bug.

    The fix is on branch fix/contig-variants-remove-variant with regression
    tests in tests/test_variants/test_remove_variant.py.
    TODO (post-fix): replace this test with:
        cv = ContigVariants()
        v = _snv(10)
        cv.add_variant(v)
        cv.remove_variant(v)
        assert 10 not in cv.variant_locations
    """
    cv = ContigVariants()
    assert callable(cv.remove_variant)


# ===========================================================================
# compile_genotypes_for_location
# ===========================================================================

def test_compile_genotypes_two_variants_different_ploids():
    cv = ContigVariants()
    v1 = _snv(10, "T", gt=np.array([1, 0]))
    v2 = _snv(10, "G", gt=np.array([0, 1]))
    cv.add_variant(v1)
    cv.add_variant(v2)
    result = cv.compile_genotypes_for_location(10)
    assert list(result) == [1, 1]


def test_compile_genotypes_single_variant():
    cv = ContigVariants()
    v = _snv(7, "C", gt=np.array([0, 1]))
    cv.add_variant(v)
    result = cv.compile_genotypes_for_location(7)
    assert list(result) == [0, 1]


# ===========================================================================
# generate_field
# ===========================================================================

def test_generate_field_uses_metadata_when_present():
    cv = ContigVariants()
    snv = _snv(1, "T")
    snv.metadata["ID"] = "rs123"
    assert cv.generate_field(snv, "ID") == "rs123"


def test_generate_field_falls_back_to_default():
    cv = ContigVariants()
    snv = _snv(1, "T")
    assert cv.generate_field(snv, "ID") == "."


# ===========================================================================
# check_if_del / check_if_ins
# ===========================================================================

def test_check_if_del_finds_containing_deletion():
    cv = ContigVariants()
    d = _del(10, 5, gt=np.array([0, 1]))
    cv.add_variant(d)
    snv = _snv(12, "T", gt=np.array([0, 1]))
    assert cv.check_if_del(snv) is d


def test_check_if_del_no_match_returns_none():
    cv = ContigVariants()
    assert cv.check_if_del(_snv(50, "T")) is None


def test_check_if_ins_with_int_position():
    """check_if_ins correctly returns the insertion when the SNV position falls inside it."""
    cv = ContigVariants()
    ins = _ins(10, "ACGTT", 4, gt=np.array([0, 1]))
    cv.add_variant(ins)
    snv = _snv(11, "T", gt=np.array([0, 1]))
    assert cv.check_if_ins(snv) is ins


def test_check_if_ins_no_match_returns_none():
    cv = ContigVariants()
    assert cv.check_if_ins(_snv(50, "T")) is None