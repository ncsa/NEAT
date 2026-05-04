"""
Tests for neat/read_simulator/utils/read.py
"""
import io

import numpy as np
import pytest
from Bio.Seq import Seq

from neat.models import SequencingErrorModel, TraditionalQualityModel
from neat.read_simulator.utils.read import Read
from neat.variants import SingleNucleotideVariant
from neat.variants.deletion import Deletion
from neat.variants.insertion import Insertion

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_READ_LEN = 100
_REF = "ACGT" * (_READ_LEN // 4)   # 100-base reference with no Ns
_PADDED_REF = _REF + "ACGT" * 5    # 120-base ref with 20 bases of padding


def _make_read(
    position=0,
    end_point=None,
    reference=_REF,
    ref_id="chr1",
    ref_id_index=0,
    padding=20,
    read_len=_READ_LEN,
    is_reverse=False,
    is_paired=False,
    raw_read=None,
):
    if end_point is None:
        end_point = position + read_len
    if raw_read is None:
        raw_read = (position, end_point, position + 150, end_point + 150)
    return Read(
        name="test_read",
        raw_read=raw_read,
        reference_segment=Seq(reference),
        reference_id=ref_id,
        ref_id_index=ref_id_index,
        position=position,
        end_point=end_point,
        padding=padding,
        run_read_len=read_len,
        is_reverse=is_reverse,
        is_paired=is_paired,
    )


def _make_rng(seed=0):
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# __repr__ and __str__
# ---------------------------------------------------------------------------

def test_repr():
    r = _make_read(position=10, end_point=110)
    assert repr(r) == "chr1: 10-110"


def test_str():
    r = _make_read(position=10, end_point=110)
    assert str(r) == "chr1: 10-110"


# ---------------------------------------------------------------------------
# Comparison operators — same reference_id
# ---------------------------------------------------------------------------

def test_gt_same_chrom_true():
    r1 = _make_read(position=200)
    r2 = _make_read(position=100)
    assert r1 > r2

def test_gt_same_chrom_false():
    r1 = _make_read(position=50)
    r2 = _make_read(position=100)
    assert not (r1 > r2)

def test_ge_same_chrom():
    r1 = _make_read(position=100)
    r2 = _make_read(position=100)
    assert r1 >= r2
    r3 = _make_read(position=50)
    assert not (r3 >= r2)

def test_lt_same_chrom_true():
    r1 = _make_read(position=50)
    r2 = _make_read(position=100)
    assert r1 < r2

def test_lt_same_chrom_false():
    r1 = _make_read(position=200)
    r2 = _make_read(position=100)
    assert not (r1 < r2)

def test_le_same_chrom():
    r1 = _make_read(position=100)
    r2 = _make_read(position=100)
    assert r1 <= r2
    r3 = _make_read(position=200)
    assert not (r3 <= r2)

def test_ne_same_chrom_different_position():
    r1 = _make_read(position=100)
    r2 = _make_read(position=200)
    assert r1 != r2

def test_eq_same_chrom_same_position():
    r1 = _make_read(position=100, end_point=200)
    r2 = _make_read(position=100, end_point=200)
    assert r1 == r2

def test_eq_same_chrom_different_end():
    r1 = _make_read(position=100, end_point=200)
    r2 = _make_read(position=100, end_point=210)
    assert r1 != r2

# ---------------------------------------------------------------------------
# Comparison operators — different reference_id
# ---------------------------------------------------------------------------

def test_gt_different_chrom_returns_false():
    r1 = _make_read(position=500, ref_id="chr1")
    r2 = _make_read(position=100, ref_id="chr2")
    assert not (r1 > r2)

def test_ge_different_chrom_returns_false():
    r1 = _make_read(position=500, ref_id="chr1")
    r2 = _make_read(position=100, ref_id="chr2")
    assert not (r1 >= r2)

def test_lt_different_chrom_returns_false():
    r1 = _make_read(position=100, ref_id="chr1")
    r2 = _make_read(position=500, ref_id="chr2")
    assert not (r1 < r2)

def test_le_different_chrom_returns_false():
    r1 = _make_read(position=100, ref_id="chr1")
    r2 = _make_read(position=500, ref_id="chr2")
    assert not (r1 <= r2)

def test_ne_different_chrom_returns_true():
    r1 = _make_read(position=100, ref_id="chr1")
    r2 = _make_read(position=100, ref_id="chr2")
    assert r1 != r2

def test_eq_different_chrom_returns_false():
    r1 = _make_read(position=100, ref_id="chr1")
    r2 = _make_read(position=100, ref_id="chr2")
    assert not (r1 == r2)


# ---------------------------------------------------------------------------
# __len__
# ---------------------------------------------------------------------------

def test_len():
    r = _make_read(read_len=150)
    assert len(r) == 150


# ---------------------------------------------------------------------------
# contains
# ---------------------------------------------------------------------------

def test_contains_inside():
    r = _make_read(position=100, end_point=200)
    assert r.contains(150)
    assert r.contains(100)       # at start (inclusive)
    assert r.contains(199)       # at end - 1 (inclusive)

def test_contains_outside():
    r = _make_read(position=100, end_point=200)
    assert not r.contains(99)
    assert not r.contains(200)   # end_point is exclusive


# ---------------------------------------------------------------------------
# update_quality_array
# ---------------------------------------------------------------------------

def _read_with_quality(length=_READ_LEN):
    r = _make_read()
    r.quality_array = np.array([30] * length, dtype=float)
    return r


def test_update_quality_array_mutation():
    r = _read_with_quality()
    r.update_quality_array(1, Seq("G"), 10, "mutation", [30], quality_score=20)
    assert r.quality_array[10] == 20
    assert len(r.quality_array) == _READ_LEN


def test_update_quality_array_error_snp():
    r = _read_with_quality()
    r.update_quality_array(1, Seq("T"), 10, "error", [30, 2])
    assert r.quality_array[10] == 2   # min of quality_scores
    assert len(r.quality_array) == _READ_LEN


def test_update_quality_array_error_insertion():
    """Insertion: len(alternate) > 1 — quality array grows by len(alt)-1 - ref_length."""
    r = _read_with_quality()
    original_len = len(r.quality_array)
    # ref_length=1, alt="ATG" (len=3) → new_quality_scores has 2 entries, replaces 1 → net +1
    r.update_quality_array(1, Seq("ATG"), 10, "error", [30, 5])
    assert len(r.quality_array) == original_len + 1
    assert r.quality_array[10] == 5
    assert r.quality_array[11] == 5


def test_update_quality_array_error_deletion():
    """Deletion: ref_length > 1 and len(alternate) == 1 — quality scores removed."""
    r = _read_with_quality()
    original_len = len(r.quality_array)
    r.update_quality_array(3, Seq("A"), 10, "error", [30, 5])
    assert len(r.quality_array) == original_len - 3


# ---------------------------------------------------------------------------
# apply_mutations
# ---------------------------------------------------------------------------

def _read_for_mutations(sequence=None, padding=20):
    r = _make_read(position=0, end_point=_READ_LEN, padding=padding)
    seq = Seq(sequence or _REF)
    r.read_sequence = seq
    r.quality_array = np.array([30] * _READ_LEN, dtype=float)
    return r


def test_apply_mutations_snv():
    r = _read_for_mutations()
    snv = SingleNucleotideVariant(
        position1=50, alt=Seq("T"), genotype=np.array([1, 1]), qual_score=30
    )
    r.mutations = {50: [snv]}
    r.apply_mutations([30], _make_rng())
    assert r.read_sequence[50] == "T"


def test_apply_mutations_insertion():
    r = _read_for_mutations()
    ins = Insertion(
        position1=50, length=2, alt=Seq("AAA"),
        genotype=np.array([1, 1]), qual_score=30
    )
    r.mutations = {50: [ins]}
    original_base = str(r.read_sequence[50])
    r.apply_mutations([30], _make_rng())
    # Insertion replaces one base with the alt sequence
    assert str(r.read_sequence[50]) != original_base or len(r.read_sequence) >= _READ_LEN


def test_apply_mutations_deletion_sufficient_padding():
    r = _read_for_mutations(padding=10)
    deletion = Deletion(
        position1=50, length=3, genotype=np.array([1, 1]), qual_score=30
    )
    r.mutations = {50: [deletion]}
    r.apply_mutations([30], _make_rng())
    assert r.padding == 7   # 10 - 3


def test_apply_mutations_deletion_insufficient_padding():
    r = _read_for_mutations(padding=0)
    deletion = Deletion(
        position1=50, length=3, genotype=np.array([1, 1]), qual_score=30
    )
    r.mutations = {50: [deletion]}
    original_seq = str(r.read_sequence)
    r.apply_mutations([30], _make_rng())
    # Deletion is skipped; sequence unchanged and padding set to 0
    assert str(r.read_sequence) == original_seq
    assert r.padding == 0


def test_apply_mutations_genotype_zero_skips():
    """genotype=[0,0] means not mutated — sequence should be unchanged."""
    r = _read_for_mutations()
    snv = SingleNucleotideVariant(
        position1=50, alt=Seq("T"), genotype=np.array([0, 0]), qual_score=30
    )
    r.mutations = {50: [snv]}
    original_seq = str(r.read_sequence)
    r.apply_mutations([30], _make_rng())
    assert str(r.read_sequence) == original_seq


# ---------------------------------------------------------------------------
# calculate_flags
# ---------------------------------------------------------------------------

def test_calculate_flags_single_ended_run():
    r = _make_read(is_paired=False, is_reverse=False)
    assert r.calculate_flags(paired_ended_run=False) == 0


def test_calculate_flags_paired_forward_proper_pair():
    # paired_ended_run, is_paired, not reverse → 1 + 2 + 32 + 64 = 99
    r = _make_read(is_paired=True, is_reverse=False)
    assert r.calculate_flags(paired_ended_run=True) == 99


def test_calculate_flags_paired_reverse_proper_pair():
    # paired_ended_run, is_paired, is_reverse → 1 + 2 + 16 + 128 = 147
    r = _make_read(is_paired=True, is_reverse=True)
    assert r.calculate_flags(paired_ended_run=True) == 147


def test_calculate_flags_paired_run_mate_unmapped():
    # paired_ended_run, not is_paired, not reverse → 1 + 8 = 9
    r = _make_read(is_paired=False, is_reverse=False)
    assert r.calculate_flags(paired_ended_run=True) == 9


def test_calculate_flags_paired_run_mate_unmapped_reverse():
    # paired_ended_run, not is_paired, is_reverse → 1 + 8 + 16 = 25
    r = _make_read(is_paired=False, is_reverse=True)
    assert r.calculate_flags(paired_ended_run=True) == 25


# ---------------------------------------------------------------------------
# get_mpos
# ---------------------------------------------------------------------------

def test_get_mpos_not_paired():
    r = _make_read(is_paired=False, raw_read=(0, 100, 200, 300))
    assert r.get_mpos() == 0


def test_get_mpos_paired_forward():
    r = _make_read(is_paired=True, is_reverse=False, raw_read=(0, 100, 200, 300))
    assert r.get_mpos() == 200   # raw_read[2]


def test_get_mpos_paired_reverse():
    r = _make_read(is_paired=True, is_reverse=True, raw_read=(0, 100, 200, 300))
    assert r.get_mpos() == 0     # raw_read[0]


# ---------------------------------------------------------------------------
# get_tlen
# ---------------------------------------------------------------------------

def test_get_tlen_not_paired():
    r = _make_read(is_paired=False, raw_read=(0, 100, 200, 300))
    assert r.get_tlen() == 0


def test_get_tlen_paired_forward():
    # length = raw_read[3] - raw_read[0] + 1 = 300 - 0 + 1 = 301
    r = _make_read(is_paired=True, is_reverse=False, raw_read=(0, 100, 200, 300))
    assert r.get_tlen() == 301


def test_get_tlen_paired_reverse():
    # same length, but negative
    r = _make_read(is_paired=True, is_reverse=True, raw_read=(0, 100, 200, 300))
    assert r.get_tlen() == -301


# ---------------------------------------------------------------------------
# convert_masking
# ---------------------------------------------------------------------------

def test_convert_masking_no_ns():
    """A reference with no Ns should be unchanged."""
    r = _make_read(reference=_REF)
    r.quality_array = np.array([30] * len(_REF), dtype=float)
    qual_model = TraditionalQualityModel()
    r.convert_masking(qual_model)
    assert "N" not in str(r.reference_segment)


def test_convert_masking_replaces_ns():
    """Ns in the reference should be replaced with TTAGGG repeat bases."""
    ref_with_n = "ACGT" * 10 + "NNNN" + "ACGT" * 15
    r = _make_read(reference=ref_with_n)
    r.quality_array = np.array([30] * len(ref_with_n), dtype=float)
    qual_model = TraditionalQualityModel()
    r.convert_masking(qual_model)
    assert "N" not in str(r.reference_segment)
    # Quality at masked positions should be set to min quality
    bad_score = min(qual_model.quality_scores)
    assert all(r.quality_array[40:44] == bad_score)


# ---------------------------------------------------------------------------
# finalize_read_and_write — produce_fastq=True
# ---------------------------------------------------------------------------

def test_finalize_read_and_write_returns_error_count():
    """Return value equals the number of errors actually applied to the read."""
    r = _make_read(reference=_PADDED_REF, padding=20)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    rng = _make_rng()

    error_count = r.finalize_read_and_write(err_model, qual_model, None, 33, False, 3, rng)

    assert isinstance(error_count, int)
    assert error_count == len(r.errors)


def test_finalize_read_and_write_writes_fastq():
    r = _make_read(reference=_PADDED_REF, padding=20)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    rng = _make_rng()
    handle = io.StringIO()

    r.finalize_read_and_write(err_model, qual_model, handle, 33, True, 3, rng)

    output = handle.getvalue()
    assert output.startswith("@test_read")
    lines = output.strip().split("\n")
    assert len(lines) == 4
    assert lines[2] == "+"
    assert len(lines[1]) == _READ_LEN
    assert len(lines[3]) == _READ_LEN


def test_finalize_read_and_write_reverse_complement():
    r = _make_read(reference=_PADDED_REF, padding=20, is_reverse=True)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    rng = _make_rng()

    r.finalize_read_and_write(err_model, qual_model, None, 33, False, 3, rng)

    assert len(r.read_sequence) == _READ_LEN


def test_finalize_sets_mapping_quality():
    r = _make_read(reference=_PADDED_REF, padding=20)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    rng = _make_rng()

    r.finalize_read_and_write(err_model, qual_model, None, 33, False, 3, rng)

    assert r.mapping_quality == 70


# ---------------------------------------------------------------------------
# make_cigar
# ---------------------------------------------------------------------------

def test_make_cigar_all_match():
    """A read identical to its reference should produce an all-M cigar."""
    r = _make_read(reference=_PADDED_REF, padding=20)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    rng = _make_rng(seed=0)
    r.finalize_read_and_write(err_model, qual_model, None, 33, False, 3, rng)
    cigar = r.make_cigar()
    assert cigar.endswith("M")
    assert "I" not in cigar or "D" not in cigar  # no complex indels for a clean read


def test_make_cigar_reverse_strand():
    """make_cigar on a reverse read should return a valid cigar string."""
    r = _make_read(reference=_PADDED_REF, padding=20, is_reverse=True)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    rng = _make_rng(seed=0)
    r.finalize_read_and_write(err_model, qual_model, None, 33, False, 3, rng)
    cigar = r.make_cigar()
    assert isinstance(cigar, str)
    assert len(cigar) > 0


# ---------------------------------------------------------------------------
# segment_start — indel position in reverse reads (PR 276 regression)
#
# For read2 (reverse), the reference segment starts `padding` bases BEFORE
# self.position to allow room for deletions before reverse-complementing.
# Variant positions must be offset from segment_start, not position.
# ---------------------------------------------------------------------------

# Build a 120-base segment (read_len=100 + padding=20) for reverse-read tests.
_SEG_LEN = _READ_LEN + 20
_REV_SEG = "ACGT" * (_SEG_LEN // 4 + 1)
_REV_SEG = _REV_SEG[:_SEG_LEN]  # exactly 120 bases

_SEGMENT_START = 80   # reference coord where segment begins
_READ2_POS    = 100   # reference coord of read2's nominal start (= segment_start + padding)
_PADDING      = _READ2_POS - _SEGMENT_START  # 20


def _make_read2(segment=_REV_SEG, segment_start=_SEGMENT_START, position=_READ2_POS,
                padding=_PADDING, read_len=_READ_LEN):
    """Return a reverse read mimicking a real read2 where segment_start < position."""
    r = Read(
        name="test_read2",
        raw_read=(segment_start, segment_start + len(segment), position + 150, position + 150 + read_len),
        reference_segment=Seq(segment),
        reference_id="chr1",
        ref_id_index=0,
        position=position,
        end_point=position + read_len,
        padding=padding,
        run_read_len=read_len,
        segment_start=segment_start,
        is_reverse=True,
        is_paired=True,
    )
    r.read_sequence = Seq(segment)
    r.quality_array = np.array([30] * len(segment), dtype=float)
    return r


def test_segment_start_defaults_to_position():
    """When segment_start is omitted, it falls back to position."""
    r = _make_read(position=50)
    assert r.segment_start == 50


def test_segment_start_stored_independently_from_position():
    """segment_start is kept separate from position for reverse reads."""
    r = _make_read2()
    assert r.segment_start == _SEGMENT_START
    assert r.position == _READ2_POS
    assert r.segment_start < r.position


def test_apply_mutations_snv_reverse_read_uses_segment_start():
    """
    Regression for PR 276: an SNV at reference position V must land at index
    V - segment_start in the segment, not V - position.

    With segment_start=80 and position=100, a variant at ref pos 110 must go
    to segment index 30 (correct) not index 10 (old, wrong).
    """
    r = _make_read2()
    variant_ref_pos = 110
    correct_idx = variant_ref_pos - _SEGMENT_START  # 30
    wrong_idx    = variant_ref_pos - _READ2_POS      # 10

    # "ACGT"*30 — index 30 is 'G', index 10 is also 'G'; use alt='T' to detect placement
    assert str(r.read_sequence[correct_idx]) == "G"
    assert str(r.read_sequence[wrong_idx]) == "G"

    snv = SingleNucleotideVariant(position1=variant_ref_pos, alt=Seq("T"),
                                  genotype=np.array([1, 1]), qual_score=30)
    r.mutations = {variant_ref_pos: [snv]}
    r.apply_mutations([30], _make_rng())

    assert str(r.read_sequence[correct_idx]) == "T", "SNV not at V - segment_start"
    assert str(r.read_sequence[wrong_idx]) == "G",   "SNV incorrectly placed at V - position"


def test_apply_mutations_snv_reverse_read_variant_before_position():
    """
    A variant in the padded region (segment_start <= V < position) must also
    be placed correctly. The old code (V - position) would give a negative index.
    """
    r = _make_read2()
    variant_ref_pos = 85  # before position=100, inside segment starting at 80
    correct_idx = variant_ref_pos - _SEGMENT_START  # 5

    # "ACGT"*30 — index 5 is 'C'
    assert str(r.read_sequence[correct_idx]) == "C"

    snv = SingleNucleotideVariant(position1=variant_ref_pos, alt=Seq("T"),
                                  genotype=np.array([1, 1]), qual_score=30)
    r.mutations = {variant_ref_pos: [snv]}
    r.apply_mutations([30], _make_rng())

    assert str(r.read_sequence[correct_idx]) == "T"


def test_apply_mutations_insertion_reverse_read_correct_position():
    """Insertion in a reverse read is placed at V - segment_start."""
    r = _make_read2()
    variant_ref_pos = 115
    correct_idx = variant_ref_pos - _SEGMENT_START  # 35

    ins = Insertion(position1=variant_ref_pos, length=2, alt=Seq("TTT"),
                    genotype=np.array([1, 1]), qual_score=30)
    r.mutations = {variant_ref_pos: [ins]}
    r.apply_mutations([30], _make_rng())

    # The three bases starting at correct_idx should now be "TTT"
    assert str(r.read_sequence[correct_idx: correct_idx + 3]) == "TTT"


def test_apply_mutations_deletion_reverse_read_correct_position():
    """Deletion in a reverse read removes bases starting at V - segment_start."""
    r = _make_read2()
    variant_ref_pos = 105
    correct_idx = variant_ref_pos - _SEGMENT_START  # 25
    del_len = 3

    # "ACGT"*30: indices 25,26,27 = 'C','G','T'; index 28 = 'A'
    pre_at_28 = str(r.read_sequence[correct_idx + del_len])  # 'A'
    pre_len = len(r.read_sequence)

    deletion = Deletion(position1=variant_ref_pos, length=del_len,
                        genotype=np.array([1, 1]), qual_score=30)
    r.mutations = {variant_ref_pos: [deletion]}
    r.apply_mutations([30], _make_rng())

    # Sequence is shorter by del_len - 1 (one base kept at position, rest removed)
    assert len(r.read_sequence) == pre_len - (del_len - 1)
    # The base that was at correct_idx + del_len is now at correct_idx + 1
    assert str(r.read_sequence[correct_idx + 1]) == pre_at_28