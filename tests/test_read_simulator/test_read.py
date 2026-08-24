"""
Tests for neat/read_simulator/utils/read.py
"""
import io
import re

import numpy as np
import pytest
from Bio.Seq import Seq

from neat.models import SequencingErrorModel, TraditionalQualityModel, ErrorContainer
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
    genomic_len=None,
    adapter_seq="",
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
        genomic_len=genomic_len,
        adapter_seq=adapter_seq,
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


def test_convert_masking_default_keeps_literal_n():
    """Under the default 'exclude' policy, residual Ns are kept as literal N at floor quality."""
    ref_with_n = "ACGT" * 10 + "NNNN" + "ACGT" * 15
    r = _make_read(reference=ref_with_n)
    r.quality_array = np.array([30] * len(ref_with_n), dtype=float)
    qual_model = TraditionalQualityModel()
    r.convert_masking(qual_model)  # defaults to n_handling="exclude"
    # The N's survive as literal base calls (no fabricated sequence)...
    assert str(r.reference_segment)[40:44] == "NNNN"
    # ...and only those positions are dropped to the model floor quality.
    bad_score = min(qual_model.quality_scores)
    assert all(r.quality_array[40:44] == bad_score)
    assert all(r.quality_array[:40] == 30)
    assert all(r.quality_array[44:] == 30)


def test_convert_masking_telomere_legacy_fills_ttaggg():
    """The legacy 'telomere' policy still fills Ns with TTAGGG repeat bases at floor quality."""
    ref_with_n = "ACGT" * 10 + "NNNN" + "ACGT" * 15
    r = _make_read(reference=ref_with_n)
    r.quality_array = np.array([30] * len(ref_with_n), dtype=float)
    qual_model = TraditionalQualityModel()
    r.convert_masking(qual_model, n_handling="telomere")
    assert "N" not in str(r.reference_segment)
    # Replacement bases are drawn from the telomere repeat alphabet.
    assert set(str(r.reference_segment)[40:44]) <= set("TTAGGG")
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


def _finalize_seq_with_n(n_handling):
    """Run finalize on an N-containing read with no errors and return the emitted FASTQ seq."""
    # 120-base segment (100 read + 20 padding); N run sits inside the read window.
    ref = "ACGT" * 10 + "NNNN" + "ACGT" * 19
    r = _make_read(reference=ref, padding=20)
    err_model = SequencingErrorModel(read_length=_READ_LEN)
    qual_model = TraditionalQualityModel()
    handle = io.StringIO()
    # num_errors=0 so no sequencing errors perturb the bases — we observe masking alone.
    r.finalize_read_and_write(err_model, qual_model, handle, 33, True, 0, _make_rng(),
                              n_handling=n_handling)
    return handle.getvalue().strip().split("\n")[1]


def test_finalize_threads_n_handling_default_literal_n():
    """finalize_read_and_write defaults to exclude: N survives as a literal base in the FASTQ."""
    seq = _finalize_seq_with_n("exclude")
    assert "N" in seq


def test_finalize_threads_n_handling_telomere():
    """Passing n_handling='telomere' reaches convert_masking: N is filled, none left literal."""
    seq = _finalize_seq_with_n("telomere")
    assert "N" not in seq


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

# ===========================================================================
# 3' adapter readthrough
# ===========================================================================

_TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


# What generate_reads budgets a short insert for deletions. A short insert has no reference
# past the fragment to draw literal headroom from, so this is an allowance rather than a
# count of bases in hand: a deletion spends it, and the adapter tail makes up the difference.
_DELETION_HEADROOM = _READ_LEN // 5


def _finalize_adapter_read(genomic_len, adapter=_TRUSEQ_R1, is_reverse=False,
                           num_errors=0, seed=0, mutations=None,
                           padding=_DELETION_HEADROOM):
    """
    Finalize a short-insert read whose reference segment is exactly `genomic_len` long,
    built the way generate_reads builds one.
    """
    r = _make_read(
        reference=_REF[:genomic_len],
        padding=padding,
        end_point=genomic_len,
        is_reverse=is_reverse,
        genomic_len=genomic_len,
        adapter_seq=adapter,
    )
    if mutations:
        r.mutations = mutations
    r.finalize_read_and_write(
        SequencingErrorModel(read_length=_READ_LEN),
        TraditionalQualityModel(),
        None,
        33,
        False,
        num_errors,
        _make_rng(seed),
    )
    return r


def _deletion_at(location, length):
    """One homozygous deletion, in NEAT's VCF-style notation (see _deleted_bases)."""
    return {location: [Deletion(position1=location, length=length,
                                genotype=np.array([1, 1]), qual_score=30)]}


def _deleted_bases(length):
    """
    Bases a NEAT Deletion of `length` actually removes from a read.

    position1 is the base *before* the first deleted one, VCF-style, and apply_mutations keeps
    it as the alternate — so a length-3 deletion takes 2 bases out of the read.
    """
    return length - 1


def _cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDS])", cigar)]


def test_adapter_defaults_leave_read_untouched():
    """No adapter_seq means no readthrough state at all — the ordinary-read path is unchanged."""
    r = _make_read()
    assert r.adapter_length == 0
    assert r.genomic_length == r.run_read_length == _READ_LEN
    assert len(r.quality_array) == _READ_LEN
    assert r.make_cigar() == f"{_READ_LEN}M"


def test_adapter_pads_short_insert_to_full_read_length():
    """A short insert is padded back out to run_read_length, sequence and quality in step."""
    genomic_len = 60
    r = _finalize_adapter_read(genomic_len)

    assert r.adapter_length == _READ_LEN - genomic_len
    assert len(r.read_sequence) == _READ_LEN
    # Sequence and quality must agree exactly: a mismatch is a malformed FASTQ record, which
    # some aligner parsers silently truncate on rather than reject.
    assert len(r.quality_array) == _READ_LEN
    assert len(r.read_quality_string) == _READ_LEN


def test_adapter_tail_matches_adapter_sequence():
    """With no sequencing errors the 3' tail is the adapter verbatim, and the 5' part is genomic."""
    # Tail of 20 bases, shorter than the 33-base adapter, so no wraparound here.
    genomic_len = _READ_LEN - 20
    r = _finalize_adapter_read(genomic_len)

    tail = str(r.read_sequence[genomic_len:])
    assert len(tail) < len(_TRUSEQ_R1)
    assert tail == _TRUSEQ_R1[:len(tail)]
    assert str(r.read_sequence[:genomic_len]) == _REF[:genomic_len]


def test_adapter_sequence_is_sourced_cyclically():
    """A read can run past the end of the adapter, so adapter bases repeat rather than run out."""
    # Genomic part short enough that the tail is longer than the adapter itself.
    genomic_len = _READ_LEN - len(_TRUSEQ_R1) - 10
    r = _finalize_adapter_read(genomic_len)

    tail = str(r.read_sequence[genomic_len:])
    assert len(tail) > len(_TRUSEQ_R1)
    expected = (_TRUSEQ_R1 * 3)[:len(tail)]
    assert tail == expected


def test_adapter_cigar_soft_clips_tail_on_forward_read():
    """Forward read: adapter is at the end in reference-forward orientation, so S trails."""
    genomic_len = 60
    r = _finalize_adapter_read(genomic_len)
    assert r.make_cigar() == f"{genomic_len}M{_READ_LEN - genomic_len}S"


def test_adapter_cigar_soft_clips_lead_on_reverse_read():
    """
    Reverse read: write_bam_record flips SEQ back to reference-forward before writing, so the
    adapter — 3' as sequenced — must appear as a LEADING soft clip in the CIGAR.
    """
    genomic_len = 60
    r = _finalize_adapter_read(genomic_len, is_reverse=True)
    assert r.make_cigar() == f"{_READ_LEN - genomic_len}S{genomic_len}M"


def test_adapter_cigar_query_length_matches_sequence():
    """M+I+S must account for every base in the read, errors or not."""
    import re

    for is_reverse in (False, True):
        r = _finalize_adapter_read(60, is_reverse=is_reverse, num_errors=5)
        ops = re.findall(r"(\d+)([MIDS])", r.make_cigar())
        query_len = sum(int(n) for n, op in ops if op in "MIS")
        assert query_len == len(r.read_sequence) == _READ_LEN


def test_adapter_bases_take_substitution_errors_only():
    """
    Adapter bases are real base calls and pick up substitution noise, but never indels — an
    indel there would change the read length that the padding exists to guarantee.
    """
    genomic_len = 20
    # A run of identical bases makes any substitution obvious, and floor-quality scores make
    # substitutions near-certain rather than rare.
    r = _make_read(
        reference=_REF[:genomic_len],
        padding=0,
        end_point=genomic_len,
        genomic_len=genomic_len,
        adapter_seq="AAAAAAAAAA",
    )
    # Bin every drawn score down to Q2 (~63% error rate), so substitutions are near-certain.
    qual_model = TraditionalQualityModel(quality_bins=[2])
    r.finalize_read_and_write(
        SequencingErrorModel(read_length=_READ_LEN), qual_model, None, 33, False, 0, _make_rng(1),
    )

    tail = str(r.read_sequence[genomic_len:])
    assert len(r.read_sequence) == _READ_LEN
    assert len(tail) == _READ_LEN - genomic_len
    assert any(base != "A" for base in tail), "expected substitution noise at floor quality"


def test_adapter_read_reports_full_length():
    """len() stays the emitted read length, which the BAM writer uses to size the record."""
    r = _finalize_adapter_read(60)
    assert len(r) == _READ_LEN


def test_apply_errors_deletion_keeps_the_anchor_bases_quality():
    """
    An error deletion keeps its anchor base -- alt *is* that base, VCF-style -- so the quality
    array has to keep the anchor's score. Dropping it left the array one shorter than the
    sequence: a malformed FASTQ record, and a BAM record samtools refuses to index.

    Long latent. A full-length read draws its quality array over the reference segment including
    the deletion headroom and trims to genomic_length afterwards, which absorbed the missing
    score; a short insert never reached here at all, because its zero padding made
    get_sequencing_errors skip every deletion before one could be applied.
    """
    r = _make_read(reference=_REF)
    r.read_sequence = Seq(_REF)
    r.quality_array = np.arange(len(_REF), dtype=int)   # distinct, so the kept score is named
    anchor_score = int(r.quality_array[10])

    # ErrorContainer(type, location, length, ref, alt) -- ref spans the anchor plus the bases
    # removed, exactly as get_sequencing_errors builds it.
    r.errors = [ErrorContainer(Deletion, 10, 3, Seq(_REF[10:14]), Seq(_REF[10]))]
    r.apply_errors(TraditionalQualityModel())

    assert len(r.quality_array) == len(r.read_sequence)
    assert int(r.quality_array[10]) == anchor_score, "the anchor kept its own score"


# ===========================================================================
# Deletions inside a short insert
#
# The regression these guard: a short insert used to be built with padding=0, and both
# apply_mutations and get_sequencing_errors skip any deletion that padding cannot cover. Every
# deletion in an adapter-readthrough read was therefore dropped, silently, wherever it fell —
# the read came back byte-identical to the unmutated reference while the golden VCF still
# claimed the variant.
# ===========================================================================

@pytest.mark.parametrize("is_reverse", [False, True])
def test_short_insert_deletion_is_applied(is_reverse):
    """The deletion survives, and the adapter tail absorbs what it removed."""
    genomic_len, del_len = 60, 4
    lost = _deleted_bases(del_len)
    r = _finalize_adapter_read(
        genomic_len, is_reverse=is_reverse, mutations=_deletion_at(20, del_len),
    )

    # Present in the CIGAR, at its real length.
    assert (lost, "D") in _cigar_ops(r.make_cigar())
    # The sequencer ran its full set of cycles either way, so the read is still full length...
    assert len(r.read_sequence) == _READ_LEN
    # ...with the adapter tail grown by exactly what the deletion took out.
    assert r.genomic_length == genomic_len - lost
    assert r.adapter_length == _READ_LEN - genomic_len + lost
    # A FASTQ record whose sequence and quality disagree is malformed.
    assert len(r.quality_array) == len(r.read_sequence)
    assert len(r.read_quality_string) == len(r.read_sequence)


@pytest.mark.parametrize("is_reverse", [False, True])
def test_short_insert_deletion_cigar_accounts_for_every_base(is_reverse):
    """M+I+S must still cover the whole read once a deletion has shortened its genomic part."""
    r = _finalize_adapter_read(
        60, is_reverse=is_reverse, mutations=_deletion_at(20, 4),
    )
    ops = _cigar_ops(r.make_cigar())
    assert sum(n for n, op in ops if op in "MIS") == len(r.read_sequence) == _READ_LEN
    # The soft clip stays on the strand-correct end (see _add_adapter_soft_clip).
    assert (ops[0][1] == "S") if is_reverse else (ops[-1][1] == "S")


def test_short_insert_deletion_without_adapter_shortens_the_read():
    """
    keep_short_fragments with no adapter: nothing backfills, so the read is emitted shorter.
    That is what a shorter sequenced molecule looks like, and every length that describes the
    read has to agree on it — run_read_length is what the BAM writer sizes the record by.
    """
    genomic_len, del_len = 60, 4
    lost = _deleted_bases(del_len)
    r = _make_read(
        reference=_REF[:genomic_len], padding=_DELETION_HEADROOM, end_point=genomic_len,
        read_len=genomic_len, genomic_len=genomic_len, adapter_seq="",
    )
    r.mutations = _deletion_at(20, del_len)
    r.finalize_read_and_write(
        SequencingErrorModel(read_length=_READ_LEN), TraditionalQualityModel(),
        None, 33, False, 0, _make_rng(0),
    )

    assert r.adapter_length == 0
    assert "S" not in r.make_cigar()
    assert len(r.read_sequence) == genomic_len - lost
    assert r.genomic_length == r.run_read_length == len(r) == genomic_len - lost
    assert len(r.quality_array) == len(r.read_sequence)


def test_full_length_read_is_untouched_by_the_resync():
    """
    An ordinary read makes a deletion up from its padding — the reference just past the window —
    and still emits exactly read_len bases. None of the short-insert bookkeeping may fire.
    """
    r = _make_read(reference=_PADDED_REF, padding=20)
    r.mutations = _deletion_at(20, 4)
    r.finalize_read_and_write(
        SequencingErrorModel(read_length=_READ_LEN), TraditionalQualityModel(),
        None, 33, False, 0, _make_rng(0),
    )

    assert r.genomic_length == r.run_read_length == _READ_LEN
    assert r.adapter_length == 0
    assert len(r.read_sequence) == _READ_LEN


# --- Boundary cases -------------------------------------------------------

def test_short_insert_deletion_at_first_genomic_base():
    """A deletion anchored on a forward read's first base is representable."""
    genomic_len, del_len = 60, 4
    lost = _deleted_bases(del_len)
    r = _finalize_adapter_read(genomic_len, mutations=_deletion_at(0, del_len))

    ops = _cigar_ops(r.make_cigar())
    assert (lost, "D") in ops
    assert r.genomic_length == genomic_len - lost
    assert sum(n for n, op in ops if op in "MIS") == len(r.read_sequence) == _READ_LEN
    # SAM cannot open an alignment on a deletion, so the CIGAR must not start with one.
    assert ops[0][1] != "D"


def test_short_insert_deletion_at_a_reverse_reads_three_prime_edge():
    """
    The same deletion on a reverse read sits at the read's 3' edge, which is the CIGAR's
    *leading* edge once SEQ is flipped to reference-forward. Known limitation: the alignment
    absorbs it rather than emitting a leading D, so the CIGAR comes back all-M and POS slides
    by the deleted length, misplacing the single anchor base.

    The read itself is still right — the deletion really was applied, the adapter grew to match,
    and the record is well formed. Only the annotation of that one edge base is approximate.

    Pricing the aligner's trailing query gaps would recover the D here, but it costs far more
    than it buys: it makes short inserts carrying an *insertion* misplace 5-10x more often,
    because there the trailing template genuinely is sequence the read never reached. Measured
    over 2400-read runs, adapter+insertion misplacement went 0.08-0.17% -> 0.96-1.62%.
    """
    genomic_len, del_len = 60, 4
    lost = _deleted_bases(del_len)
    r = _finalize_adapter_read(genomic_len, is_reverse=True, mutations=_deletion_at(0, del_len))

    # Applied to the read, whatever the CIGAR says about it.
    assert r.genomic_length == genomic_len - lost
    assert r.adapter_length == _READ_LEN - genomic_len + lost
    assert len(r.read_sequence) == _READ_LEN
    assert len(r.quality_array) == len(r.read_sequence)

    cigar, reference_start = r.make_alignment()
    ops = _cigar_ops(cigar)
    assert sum(n for n, op in ops if op in "MIS") == _READ_LEN
    assert ops[0][1] == "S"          # adapter still clipped on the strand-correct end
    assert ops[1][1] != "D"          # and the alignment never opens on a deletion
    # POS absorbs the shift instead, and stays inside the read's own window.
    assert r.position <= reference_start <= r.position + lost


@pytest.mark.parametrize("is_reverse", [False, True])
def test_short_insert_deletion_at_final_genomic_base(is_reverse):
    """
    A deletion anchored on the last base of the fragment. Its anchor is kept and the bases it
    would remove lie past the molecule, so nothing is removed — but it must not corrupt the
    record on the way through.
    """
    genomic_len = 60
    r = _finalize_adapter_read(
        genomic_len, is_reverse=is_reverse,
        mutations=_deletion_at(genomic_len - 1, 4),
    )

    assert r.genomic_length == genomic_len
    assert len(r.read_sequence) == _READ_LEN
    assert len(r.quality_array) == len(r.read_sequence)
    ops = _cigar_ops(r.make_cigar())
    assert sum(n for n, op in ops if op in "MIS") == _READ_LEN


@pytest.mark.parametrize("is_reverse", [False, True])
def test_short_insert_deletion_extending_past_the_fragment(is_reverse):
    """
    A deletion whose span runs off the 3' end of the fragment. Only the part inside the molecule
    can be sequenced away, so it is applied up to the edge and no further — and the read stays
    well formed, which is the property that matters for the golden BAM.
    """
    genomic_len, del_len = 60, 10
    overhang_anchor = genomic_len - 4          # 6 of the 9 deleted bases lie past the fragment
    r = _finalize_adapter_read(
        genomic_len, is_reverse=is_reverse, mutations=_deletion_at(overhang_anchor, del_len),
    )

    # Clamped at the fragment edge: at most the bases that were actually there.
    assert genomic_len - _deleted_bases(del_len) <= r.genomic_length < genomic_len
    assert r.genomic_length == overhang_anchor + 1
    assert len(r.read_sequence) == _READ_LEN
    assert len(r.quality_array) == len(r.read_sequence)
    ops = _cigar_ops(r.make_cigar())
    assert sum(n for n, op in ops if op in "MIS") == _READ_LEN


def test_deletion_with_no_headroom_left_is_logged(caplog):
    """
    The skip is a last resort — it drops ground truth the golden VCF still carries — so it may
    not be silent. Two deletions, the second of which the budget cannot cover.
    """
    genomic_len = 60
    mutations = _deletion_at(10, _DELETION_HEADROOM)
    mutations.update(_deletion_at(30, 8))
    with caplog.at_level("DEBUG", logger="neat.read_simulator.utils.read"):
        _finalize_adapter_read(genomic_len, mutations=mutations)

    assert any("Skipped a" in message for message in caplog.messages)


# ---------------------------------------------------------------------------
# CIGAR construction via alignment (issue 326)
#
# A non-repetitive reference: _REF is "ACGT" repeated, which an aligner can
# match at several offsets, so these tests supply their own template.
# ---------------------------------------------------------------------------

_UNIQUE_REF = (
    "TTGACCATGGCAGTTCAAGGCTATCCGAATTCACGGTACCTAGGCATTAGCCGGATCAATGCC"
    "AAGGTTCCAATGGCTTAAGCCTGATCAGGTTACCGGAATTCCGGTTAACCGGATTACGGCATA"
)


def _alignment_read(read_sequence, reference, position=0, end_point=None,
                    is_reverse=False, genomic_len=None, segment_start=None):
    """A read whose CIGAR must come from the alignment path, with its sequence supplied."""
    genomic_len = genomic_len if genomic_len is not None else len(read_sequence)
    if end_point is None:
        end_point = position + genomic_len
    r = Read(
        name="align_read",
        raw_read=(position, end_point, position + 150, end_point + 150),
        reference_segment=Seq(reference),
        reference_id="chr1",
        ref_id_index=0,
        position=position,
        end_point=end_point,
        padding=len(reference) - genomic_len,
        run_read_len=genomic_len,
        segment_start=segment_start if segment_start is not None else position,
        is_reverse=is_reverse,
        genomic_len=genomic_len,
    )
    r.read_sequence = Seq(read_sequence)
    # A mutation indel on the read is what routes make_cigar to the alignment path.
    r.mutations = {position: [Insertion(position, 2, "AT", np.array([1, 1]))]}
    return r


def _ops(cigar):
    import re
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([MIDS])", cigar)]


def test_reference_span_counts_only_reference_consuming_ops():
    assert Read.reference_span("100M") == 100
    assert Read.reference_span("10M5I85M") == 95
    assert Read.reference_span("10M5D85M") == 100
    assert Read.reference_span("10S90M") == 90
    assert Read.reference_span("50M600I") == 50


def test_long_insertion_survives_into_the_cigar():
    """An insertion longer than 4 bp used to be unrepresentable: the op list was fixed at the
    read length and an insertion was recorded by overwriting one of its entries."""
    anchor = _UNIQUE_REF[:40]
    inserted = "GGGGTTTTGGGGTTTTGGGGTTTTGGGGTTTTGGGGTTTTGGGGTTTTGGGGTTTTGGGG"
    read = _alignment_read(anchor + inserted, _UNIQUE_REF)
    cigar, _ = read.make_alignment()
    ops = _ops(cigar)
    # Not the full 60: the tail of a repetitive insert can find a match in the template. The point
    # is the order of magnitude — the old op list could not carry an insertion past 4 bp at all.
    assert max((n for op, n in ops if op == "I"), default=0) >= 40
    # every base of the read is still accounted for
    assert sum(n for op, n in ops if op in "MIS") == len(read.read_sequence)


def test_insertion_longer_than_the_read_is_bounded_by_it():
    """A read made entirely of inserted sequence cannot claim more query bases than it has."""
    inserted = "GGGGTTTT" * 12
    read = _alignment_read(inserted[:64], _UNIQUE_REF)
    cigar, _ = read.make_alignment()
    ops = _ops(cigar)
    assert sum(n for op, n in ops if op in "MIS") == 64


def test_cigar_query_length_always_matches_the_read():
    """M+I+S must cover the read whichever path built the cigar."""
    for cut in (10, 40, 80):
        read = _alignment_read(_UNIQUE_REF[:cut] + "GGGGTTTT" * 4, _UNIQUE_REF)
        ops = _ops(read.make_cigar())
        assert sum(n for op, n in ops if op in "MIS") == len(read.read_sequence)


def test_reverse_read_position_is_anchored_on_its_right_edge():
    """A reverse read's headroom sits before its window, so a deletion moves where the alignment
    starts. POS has to come from end_point minus the reference the cigar covers, not from
    self.position, which left every gapped reverse read shifted by its net indel length."""
    headroom = 10
    # The 101 reference bases the read covers: one base of headroom plus its window. Dropping an
    # interior base leaves a 100 bp read spanning 101 bp of reference, so it reaches one base
    # further left than self.position.
    covered = _UNIQUE_REF[headroom - 1:headroom + 100]
    mutated = covered[:51] + covered[52:]
    read = _alignment_read(
        str(Seq(mutated).reverse_complement()),
        _UNIQUE_REF[:headroom + 100],
        position=headroom,
        end_point=headroom + 100,
        is_reverse=True,
        segment_start=0,
    )
    cigar, reference_start = read.make_alignment()
    assert "D" in cigar
    # A reverse read's alignment ends where its window ends, and the deletion pushes its start left.
    assert reference_start + Read.reference_span(cigar) == read.end_point
    assert reference_start == read.position - 1


def test_forward_read_position_is_anchored_on_its_left_edge():
    read = _alignment_read(_UNIQUE_REF[:40] + "GGGGTTTT" * 4, _UNIQUE_REF)
    _, reference_start = read.make_alignment()
    assert reference_start == read.position


def test_ungapped_read_reports_its_own_position():
    """The fast path must not disturb either anchor."""
    for is_reverse in (False, True):
        r = _make_read(reference=_PADDED_REF, padding=20, is_reverse=is_reverse)
        r.read_sequence = Seq(_REF)
        cigar, reference_start = r.make_alignment()
        assert cigar == f"{_READ_LEN}M"
        assert reference_start == r.position
