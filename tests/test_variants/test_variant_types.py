"""
Tests for variant type classes: Deletion, Insertion, SingleNucleotideVariant, UnknownVariant.
"""

import numpy as np
import pytest
from neat.variants import Deletion, Insertion, SingleNucleotideVariant
from neat.variants.unknown_variant import UnknownVariant


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

GENOTYPE_HET = np.array([0, 1])
GENOTYPE_HOM = np.array([1, 1])
GENOTYPE_REF = np.array([0, 0])


# ===========================================================================
# Deletion tests
# ===========================================================================

class TestDeletion:
    """Tests for Deletion."""

    def _make(self, position1=100, length=5, genotype=None, qual_score=30, is_input=False, **kwargs):
        if genotype is None:
            genotype = GENOTYPE_HET.copy()
        return Deletion(position1, length, genotype, qual_score, is_input, **kwargs)

    # --- __repr__ -----------------------------------------------------------

    def test_repr(self):
        d = self._make(position1=10, length=3)
        assert repr(d) == "Deletion(10, 3)"

    # --- get_alt ------------------------------------------------------------

    def test_get_alt_returns_empty_string(self):
        d = self._make()
        assert d.get_alt() == ""

    # --- contains -----------------------------------------------------------

    def test_contains_true_at_start(self):
        d = self._make(position1=100, length=5, genotype=GENOTYPE_HET.copy())
        other = type("Fake", (), {"position1": 100, "genotype": GENOTYPE_HET.copy()})()
        assert d.contains(other) is True

    def test_contains_true_inside(self):
        d = self._make(position1=100, length=5, genotype=GENOTYPE_HET.copy())
        other = type("Fake", (), {"position1": 103, "genotype": GENOTYPE_HET.copy()})()
        assert d.contains(other) is True

    def test_contains_false_at_end(self):
        """position1 == position1 + length is exclusive upper bound."""
        d = self._make(position1=100, length=5, genotype=GENOTYPE_HET.copy())
        other = type("Fake", (), {"position1": 105, "genotype": GENOTYPE_HET.copy()})()
        assert d.contains(other) is False

    def test_contains_false_before_range(self):
        d = self._make(position1=100, length=5, genotype=GENOTYPE_HET.copy())
        other = type("Fake", (), {"position1": 99, "genotype": GENOTYPE_HET.copy()})()
        assert d.contains(other) is False

    def test_contains_false_after_range(self):
        d = self._make(position1=100, length=5, genotype=GENOTYPE_HET.copy())
        other = type("Fake", (), {"position1": 110, "genotype": GENOTYPE_HET.copy()})()
        assert d.contains(other) is False

    def test_contains_false_genotype_differs(self):
        d = self._make(position1=100, length=5, genotype=GENOTYPE_HET.copy())
        other = type("Fake", (), {"position1": 102, "genotype": GENOTYPE_HOM.copy()})()
        assert d.contains(other) is False

    # --- comparison operators -----------------------------------------------

    def test_lt_true(self):
        """position < self.position1 → __lt__ is True"""
        d = self._make(position1=100, length=5)
        assert d.__lt__(99) is True

    def test_lt_false(self):
        d = self._make(position1=100, length=5)
        assert d.__lt__(100) is False

    def test_gt_true(self):
        """position > self.position1 + self.length → __gt__ is True"""
        d = self._make(position1=100, length=5)
        assert d.__gt__(106) is True

    def test_gt_false(self):
        d = self._make(position1=100, length=5)
        assert d.__gt__(105) is False

    def test_le_true(self):
        """position <= self.position1 → __le__ is True"""
        d = self._make(position1=100, length=5)
        assert d.__le__(100) is True
        assert d.__le__(99) is True

    def test_le_false(self):
        d = self._make(position1=100, length=5)
        assert d.__le__(101) is False

    def test_ge_true(self):
        """position >= self.position1 + self.length → __ge__ is True"""
        d = self._make(position1=100, length=5)
        assert d.__ge__(105) is True
        assert d.__ge__(110) is True

    def test_ge_false(self):
        d = self._make(position1=100, length=5)
        assert d.__ge__(104) is False

    # --- __eq__ -------------------------------------------------------------

    def test_eq_matching(self):
        """Two Deletion objects with matching type/position/length."""
        d = self._make(position1=50, length=4)
        other = type("Fake", (), {
            "type": Deletion,
            "position": 50,
            "length": 4,
        })()
        assert d.__eq__(other) is True

    def test_eq_position_differs(self):
        d = self._make(position1=50, length=4)
        other = type("Fake", (), {
            "type": Deletion,
            "position": 51,
            "length": 4,
        })()
        assert d.__eq__(other) is False

    def test_eq_length_differs(self):
        d = self._make(position1=50, length=4)
        other = type("Fake", (), {
            "type": Deletion,
            "position": 50,
            "length": 5,
        })()
        assert d.__eq__(other) is False

    def test_eq_wrong_type(self):
        d = self._make(position1=50, length=4)
        other = type("Fake", (), {
            "type": Insertion,
            "position": 50,
            "length": 4,
        })()
        assert d.__eq__(other) is False

    # --- BaseVariant methods -------------------------------------------------

    def test_get_qual_score_direct(self):
        d = self._make(qual_score=42)
        assert d.get_qual_score() == 42

    def test_get_qual_score_from_metadata(self):
        d = self._make(qual_score=None, QUAL=99)
        assert d.get_qual_score() == 99

    def test_get_0_location(self):
        d = self._make(position1=77)
        assert d.get_0_location() == 77

    def test_get_1_location(self):
        d = self._make(position1=77)
        assert d.get_1_location() == 78


# ===========================================================================
# Insertion tests
# ===========================================================================

class TestInsertion:
    """Tests for Insertion."""

    def _make(self, position1=200, length=4, alt="ACGT", genotype=None,
              qual_score=30, is_input=False, **kwargs):
        if genotype is None:
            genotype = GENOTYPE_HET.copy()
        return Insertion(position1, length, alt, genotype, qual_score, is_input, **kwargs)

    # --- __repr__ -----------------------------------------------------------

    def test_repr(self):
        ins = self._make(position1=20, alt="ATCG")
        assert repr(ins) == "Insertion(20, ATCG)"

    # --- contains -----------------------------------------------------------

    def test_contains_true_at_start(self):
        ins = self._make(position1=200, length=4)
        assert ins.contains(200) is True

    def test_contains_true_inside(self):
        ins = self._make(position1=200, length=4)
        assert ins.contains(202) is True

    def test_contains_false_at_end(self):
        ins = self._make(position1=200, length=4)
        assert ins.contains(204) is False

    def test_contains_false_before(self):
        ins = self._make(position1=200, length=4)
        assert ins.contains(199) is False

    def test_contains_false_after(self):
        ins = self._make(position1=200, length=4)
        assert ins.contains(210) is False

    # --- comparison operators -----------------------------------------------

    def test_lt_true(self):
        ins = self._make(position1=200, length=4)
        assert ins.__lt__(199) is True

    def test_lt_false(self):
        ins = self._make(position1=200, length=4)
        assert ins.__lt__(200) is False

    def test_gt_true(self):
        ins = self._make(position1=200, length=4)
        assert ins.__gt__(205) is True

    def test_gt_false(self):
        ins = self._make(position1=200, length=4)
        assert ins.__gt__(204) is False

    def test_le_true(self):
        ins = self._make(position1=200, length=4)
        assert ins.__le__(200) is True
        assert ins.__le__(199) is True

    def test_le_false(self):
        ins = self._make(position1=200, length=4)
        assert ins.__le__(201) is False

    def test_ge_true(self):
        ins = self._make(position1=200, length=4)
        assert ins.__ge__(204) is True

    def test_ge_false(self):
        ins = self._make(position1=200, length=4)
        assert ins.__ge__(203) is False

    # --- __eq__ -------------------------------------------------------------

    def test_eq_matching(self):
        ins = self._make(position1=200, length=4, alt="ACGT")
        other = type("Fake", (), {
            "type": Insertion,
            "position": 200,
            "alt": "ACGT",
            "length": 4,
        })()
        assert ins.__eq__(other) is True

    def test_eq_wrong_type(self):
        ins = self._make(position1=200, length=4, alt="ACGT")
        other = type("Fake", (), {
            "type": Deletion,
            "position": 200,
            "alt": "ACGT",
            "length": 4,
        })()
        assert ins.__eq__(other) is False

    # --- BaseVariant methods via Insertion -----------------------------------

    def test_get_alt_direct(self):
        ins = self._make(alt="TTTT")
        assert ins.get_alt() == "TTTT"

    def test_get_alt_from_metadata(self):
        ins = self._make(alt=None, ALT="GGGG")
        assert ins.get_alt() == "GGGG"

    def test_get_0_location(self):
        ins = self._make(position1=333)
        assert ins.get_0_location() == 333

    def test_get_1_location(self):
        ins = self._make(position1=333)
        assert ins.get_1_location() == 334


# ===========================================================================
# SingleNucleotideVariant tests
# ===========================================================================

class TestSingleNucleotideVariant:
    """Tests for SingleNucleotideVariant."""

    def _make(self, position1=50, alt="T", genotype=None,
              qual_score=30, is_input=False, **kwargs):
        if genotype is None:
            genotype = GENOTYPE_HET.copy()
        return SingleNucleotideVariant(position1, alt, genotype, qual_score, is_input, **kwargs)

    # --- __repr__ -----------------------------------------------------------

    def test_repr(self):
        snv = self._make(position1=5, alt="G")
        assert repr(snv) == "SingleNucleotideVariant(5, G)"

    # --- comparison operators -----------------------------------------------

    def test_lt_true(self):
        snv = self._make(position1=50)
        assert snv.__lt__(49) is True

    def test_lt_false(self):
        snv = self._make(position1=50)
        assert snv.__lt__(50) is False

    def test_gt_true(self):
        snv = self._make(position1=50)
        assert snv.__gt__(51) is True

    def test_gt_false(self):
        snv = self._make(position1=50)
        assert snv.__gt__(50) is False

    def test_le_true(self):
        snv = self._make(position1=50)
        assert snv.__le__(50) is True
        assert snv.__le__(49) is True

    def test_le_false(self):
        snv = self._make(position1=50)
        assert snv.__le__(51) is False

    def test_ge_true(self):
        snv = self._make(position1=50)
        assert snv.__ge__(50) is True
        assert snv.__ge__(51) is True

    def test_ge_false(self):
        snv = self._make(position1=50)
        assert snv.__ge__(49) is False

    # --- __eq__ -------------------------------------------------------------

    def test_eq_matching(self):
        snv = self._make(position1=50, alt="T")
        other = SingleNucleotideVariant(50, "T", GENOTYPE_HET.copy(), 30, False)
        assert snv.__eq__(other) is True

    def test_eq_position_differs(self):
        snv = self._make(position1=50, alt="T")
        other = SingleNucleotideVariant(51, "T", GENOTYPE_HET.copy(), 30, False)
        assert snv.__eq__(other) is False

    def test_eq_alt_differs(self):
        snv = self._make(position1=50, alt="T")
        other = SingleNucleotideVariant(50, "A", GENOTYPE_HET.copy(), 30, False)
        assert snv.__eq__(other) is False

    def test_eq_wrong_type(self):
        snv = self._make(position1=50, alt="T")
        other = Deletion(50, 1, GENOTYPE_HET.copy(), 30, False)
        assert snv.__eq__(other) is False


# ===========================================================================
# UnknownVariant tests
# ===========================================================================

class TestUnknownVariant:
    """Tests for UnknownVariant."""

    def _make(self, position1=300, genotype=None, qual_score=20,
              is_input=True, **kwargs):
        if genotype is None:
            genotype = GENOTYPE_HET.copy()
        return UnknownVariant(position1, genotype, qual_score, is_input, **kwargs)

    # --- __repr__ -----------------------------------------------------------

    def test_repr(self):
        uv = self._make(position1=42)
        assert repr(uv) == "UnknownVariant(42)"

    # --- get_ref_len --------------------------------------------------------

    def test_get_ref_len(self):
        uv = self._make(REF="ACGT")
        assert uv.get_ref_len() == 4

    def test_get_ref_len_single(self):
        uv = self._make(REF="A")
        assert uv.get_ref_len() == 1

    # --- comparison operators -----------------------------------------------

    def test_lt_true(self):
        uv = self._make(position1=300)
        assert uv.__lt__(299) is True

    def test_lt_false(self):
        uv = self._make(position1=300)
        assert uv.__lt__(300) is False

    def test_gt_true(self):
        uv = self._make(position1=300)
        assert uv.__gt__(301) is True

    def test_gt_false(self):
        uv = self._make(position1=300)
        assert uv.__gt__(300) is False

    def test_le_true(self):
        uv = self._make(position1=300)
        assert uv.__le__(300) is True
        assert uv.__le__(299) is True

    def test_le_false(self):
        uv = self._make(position1=300)
        assert uv.__le__(301) is False

    def test_ge_true(self):
        uv = self._make(position1=300)
        assert uv.__ge__(300) is True
        assert uv.__ge__(301) is True

    def test_ge_false(self):
        uv = self._make(position1=300)
        assert uv.__ge__(299) is False

    # --- __eq__ -------------------------------------------------------------

    def test_eq_wrong_type_returns_false(self):
        """UnknownVariant.__eq__ checks other.type; a plain SNV object has no .type attr,
        which means accessing other.type raises AttributeError. We use a fake with wrong type."""
        uv = self._make(position1=300)
        other = type("Fake", (), {
            "type": Deletion,
            "position1": 300,
            "alt": None,
            "genotype": GENOTYPE_HET.copy(),
        })()
        assert uv.__eq__(other) is False

    # --- BaseVariant methods via UnknownVariant -----------------------------

    def test_get_qual_score_direct(self):
        uv = self._make(qual_score=55)
        assert uv.get_qual_score() == 55

    def test_get_qual_score_from_metadata(self):
        uv = self._make(qual_score=None, QUAL=77)
        assert uv.get_qual_score() == 77

    def test_get_0_location(self):
        uv = self._make(position1=400)
        assert uv.get_0_location() == 400

    def test_get_1_location(self):
        uv = self._make(position1=400)
        assert uv.get_1_location() == 401