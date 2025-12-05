"""
Unit tests for MutationModel and SequencingErrorModel
"""

import numpy as np
from numpy.random import default_rng
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.models.mutation_model import MutationModel
from neat.models.error_models import SequencingErrorModel, TraditionalQualityModel
from neat.variants.insertion import Insertion
from neat.variants.deletion import Deletion
from neat.variants.single_nucleotide_variant import SingleNucleotideVariant

def test_mutation_model_generate_insertion_basic():
    rng = default_rng(1)
    m = MutationModel()
    ref = Seq("A")
    ins = m.generate_insertion(location=10, ref=ref, rng=rng)
    assert isinstance(ins, Insertion)
    assert ins.length >= 1 # length stored is VCF-style (>=1)
    assert isinstance(ins.alt, (str, Seq)) # alt should be ref + inserted sequence
    assert len(str(ins.alt)) >= len(str(ref)) + 1 # inserted >= 1 base


def test_mutation_model_generate_deletion_basic():
    rng = default_rng(2)
    m = MutationModel()
    dele = m.generate_deletion(location=15, rng=rng)
    assert isinstance(dele, Deletion)
    assert dele.length >= 1 # VCF-style includes shared base
    assert hasattr(dele, "position1")


def test_mutation_model_generate_snv_trinuc():
    rng = default_rng(3)
    m = MutationModel()
    snv = m.generate_snv(trinucleotide=Seq("ACA"), reference_location=5, rng=rng)
    assert isinstance(snv, SingleNucleotideVariant)
    assert snv.alt in ["A", "C", "G", "T"]


def test_sequencing_error_model_zero_error_returns_none_or_empty():
    """
    avg_seq_error == 0 should yield no errors.
    """
    rng = default_rng(4)
    sem = SequencingErrorModel(avg_seq_error=0.0)
    ref = SeqRecord(Seq("A" * 40), id="chr1")
    quals = np.array([40] * 40, dtype=int)
    result = sem.get_sequencing_errors(
        padding=20,
        reference_segment=ref,
        quality_scores=quals,
        rng=rng,
    )
    if isinstance(result, tuple):
        introduced, pad = result
        assert introduced == []
        assert pad >= 0
    else:
        assert result == []


def test_traditional_quality_model_shapes_and_range():
    rng = default_rng(5)
    qm = TraditionalQualityModel(average_error=0.01)
    qs = qm.get_quality_scores(model_read_length=151, length=75, rng=rng)
    assert isinstance(qs, np.ndarray)
    assert len(qs) == 75
    assert qs.min() >= 1 and qs.max() <= 42


def test_sequencing_error_model_basic_snvs_only():
    rng = default_rng(6)
    sem = SequencingErrorModel(avg_seq_error=0.05)
    ref = SeqRecord(Seq("ACGT" * 20), id="chr1")
    quals = np.array([35] * 80, dtype=int)
    introduced, pad = sem.get_sequencing_errors(
        padding=40, reference_segment=ref, quality_scores=quals, rng=rng
    )
    assert isinstance(introduced, list)
    assert pad >= 0
    # If any errors occurred, they must have expected attributes
    for e in introduced:
        assert hasattr(e, "error_type")
        assert hasattr(e, "location")
        assert hasattr(e, "ref")
        assert hasattr(e, "alt")

def test_mutation_model_insertion_reproducible_with_seed():
    """Same seed and inputs should give the same insertion (length and alt)."""
    rng1 = default_rng(123)
    rng2 = default_rng(123)
    m = MutationModel()
    ref = Seq("ACGT")

    ins1 = m.generate_insertion(location=10, ref=ref, rng=rng1)
    ins2 = m.generate_insertion(location=10, ref=ref, rng=rng2)

    assert isinstance(ins1, Insertion)
    assert isinstance(ins2, Insertion)
    assert ins1.length == ins2.length
    assert str(ins1.alt) == str(ins2.alt)


def test_mutation_model_deletion_reproducible_with_seed():
    """Same seed and inputs should give the same deletion object shape."""
    rng1 = default_rng(456)
    rng2 = default_rng(456)
    m = MutationModel()

    del1 = m.generate_deletion(location=25, rng=rng1)
    del2 = m.generate_deletion(location=25, rng=rng2)

    assert isinstance(del1, Deletion)
    assert isinstance(del2, Deletion)
    assert del1.length == del2.length
    assert del1.position1 == del2.position1


def test_mutation_model_snv_does_not_keep_reference_base():
    """
    For a given trinucleotide, the generated SNV should change the central base.
    """
    rng = default_rng(7)
    m = MutationModel()
    trinuc = Seq("ACA")
    central = str(trinuc[1])

    snv = m.generate_snv(trinucleotide=trinuc, reference_location=100, rng=rng)
    assert isinstance(snv, SingleNucleotideVariant)
    assert snv.alt in ["A", "C", "G", "T"]
    assert snv.alt != central


def test_traditional_quality_model_reproducible_with_seed():
    """Quality model should be deterministic given the same RNG state."""
    rng1 = default_rng(8)
    rng2 = default_rng(8)
    qm = TraditionalQualityModel(average_error=0.01)

    qs1 = qm.get_quality_scores(model_read_length=151, length=100, rng=rng1)
    qs2 = qm.get_quality_scores(model_read_length=151, length=100, rng=rng2)

    assert np.array_equal(qs1, qs2)


def test_sequencing_error_model_reproducible_with_seed():
    """Error placement should be deterministic given the same RNG state."""
    sem = SequencingErrorModel(avg_seq_error=0.05)
    ref = SeqRecord(Seq("ACGT" * 30), id="chr1")
    quals = np.array([30] * len(ref), dtype=int)

    rng1 = default_rng(9)
    rng2 = default_rng(9)

    introduced1, pad1 = sem.get_sequencing_errors(
        padding=20, reference_segment=ref, quality_scores=quals, rng=rng1
    )
    introduced2, pad2 = sem.get_sequencing_errors(
        padding=20, reference_segment=ref, quality_scores=quals, rng=rng2
    )

    assert pad1 == pad2
    proj1 = [(e.error_type, e.location, e.ref, e.alt) for e in introduced1]
    proj2 = [(e.error_type, e.location, e.ref, e.alt) for e in introduced2]
    assert proj1 == proj2


def test_sequencing_error_model_nonzero_error_introduces_in_bounds_errors():
    """
    With a non-zero average error rate, we expect at least some errors and their
    locations must be within the reference segment.
    """
    rng = default_rng(10)
    sem = SequencingErrorModel(avg_seq_error=0.2)
    ref = SeqRecord(Seq("ACGT" * 50), id="chr1")
    quals = np.array([10] * len(ref), dtype=int)

    introduced, pad = sem.get_sequencing_errors(
        padding=20, reference_segment=ref, quality_scores=quals, rng=rng
    )

    # At least one error is expected for these settings.
    assert len(introduced) > 0
    # All error locations should be within the sequence.
    for e in introduced:
        assert 0 <= e.location < len(ref)
        assert e.ref in ["A", "C", "G", "T"]
        assert e.alt in ["A", "C", "G", "T"]
