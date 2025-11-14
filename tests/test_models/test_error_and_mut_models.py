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
    avg_seq_error == 0 should yield no errors. Some versions return just the list,
    others return a (list, padding) tuple — accept both.
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
