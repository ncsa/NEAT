import numpy as np
import pytest
from pathlib import Path
import pickle
import os
from neat.models.gc_bias_model import GCBiasModel, get_uniform_gc_model

def test_gc_bias_model_init():
    weights = [1.0] * 101
    model = GCBiasModel(weights, window_size=100)
    assert model.window_size == 100
    assert len(model.weights) == 101
    assert model.is_uniform
    assert model.max_weight == 1.0

def test_gc_bias_model_invalid_init():
    with pytest.raises(ValueError):
        GCBiasModel([1.0] * 100, 100)
    with pytest.raises(ValueError):
        GCBiasModel([1.0] * 101, 0)
    with pytest.raises(ValueError):
        GCBiasModel([-1.0] * 101, 100)

def test_gc_bias_model_get_weight():
    weights = [0.0] * 101
    weights[50] = 1.0
    model = GCBiasModel(weights, window_size=100)
    
    assert model.get_weight(0.5) == 1.0
    assert model.get_weight(0.496) == 1.0 # rounds to 0.50
    assert model.get_weight(0.0) == 0.0
    assert model.get_weight(1.0) == 0.0

def test_gc_bias_model_get_weight_for_sequence():
    weights = [0.0] * 101
    weights[50] = 1.0
    model = GCBiasModel(weights, window_size=100)
    
    # 50% GC
    assert model.get_weight_for_sequence("GCAT") == 1.0
    # 100% GC
    assert model.get_weight_for_sequence("GGCC") == weights[100]
    # 0% GC
    assert model.get_weight_for_sequence("AATT") == weights[0]
    
    # Empty/N
    assert model.get_weight_for_sequence("NNNN") == 1.0

def test_gc_bias_model_save_load(tmp_path):
    weights = np.random.random(101).tolist()
    window_size = 150
    model = GCBiasModel(weights, window_size)
    
    path = tmp_path / "gc_model.pickle"
    model.save(path)
    
    loaded_model = GCBiasModel.from_file(path)
    assert loaded_model.window_size == window_size
    np.testing.assert_allclose(loaded_model.weights, weights)

def test_uniform_gc_model():
    model = get_uniform_gc_model(200)
    assert model.window_size == 200
    assert model.is_uniform
    assert model.get_weight(0.3) == 1.0


# ---------------------------------------------------------------------------
# get_gc_bias_weights
# ---------------------------------------------------------------------------

import pysam
from neat.model_gc_bias.utils import get_gc_bias_weights


def _make_biased_bam(bam_path, read_count_at=5, read_count_gc=200):
    """Write a tiny BAM: 10 000 bp ref (5000 A + 5000 G), biased toward GC region."""
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10000, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(read_count_at):
            a = pysam.AlignedSegment()
            a.query_name = f"at_{i}"
            a.query_sequence = "A" * 50
            a.reference_id = 0
            a.reference_start = 500 + i
            a.cigar = ((0, 50),)
            a.qual = "B" * 50
            bam.write(a)
        for i in range(read_count_gc):
            a = pysam.AlignedSegment()
            a.query_name = f"gc_{i}"
            a.query_sequence = "G" * 50
            a.reference_id = 0
            a.reference_start = 6000 + i
            a.cigar = ((0, 50),)
            a.qual = "B" * 50
            bam.write(a)
    pysam.index(str(bam_path))


def test_get_gc_bias_weights_returns_101_values(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 5000 + "G" * 5000 + "\n")
    bam = tmp_path / "in.bam"
    _make_biased_bam(bam)

    weights = get_gc_bias_weights(str(bam), str(ref), window_size=100)
    assert len(weights) == 101


def test_get_gc_bias_weights_all_positive(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 5000 + "G" * 5000 + "\n")
    bam = tmp_path / "in.bam"
    _make_biased_bam(bam)

    weights = get_gc_bias_weights(str(bam), str(ref), window_size=100)
    assert all(w > 0 for w in weights), "Every weight should be positive after fill"


def test_get_gc_bias_weights_max_is_one(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 5000 + "G" * 5000 + "\n")
    bam = tmp_path / "in.bam"
    _make_biased_bam(bam)

    weights = get_gc_bias_weights(str(bam), str(ref), window_size=100)
    assert abs(max(weights) - 1.0) < 1e-9


def test_get_gc_bias_weights_gc_rich_heavier_than_at_rich(tmp_path):
    """100% GC region (bin 100) should have higher weight than 0% GC region (bin 0)."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 5000 + "G" * 5000 + "\n")
    bam = tmp_path / "in.bam"
    _make_biased_bam(bam, read_count_at=5, read_count_gc=200)

    weights = get_gc_bias_weights(str(bam), str(ref), window_size=100)
    assert weights[100] > weights[0]


# ---------------------------------------------------------------------------
# compute_genome_wide_gc_mean_weight
# ---------------------------------------------------------------------------

from neat.model_gc_bias.utils import compute_genome_wide_gc_mean_weight


def test_compute_genome_wide_gc_mean_uniform_model(tmp_path):
    """Uniform model: returns weights[0] without scanning the reference."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "ACGT" * 1000 + "\n")
    model = GCBiasModel([0.7] * 101, window_size=100)
    # No need to even index — uniform short-circuits
    assert compute_genome_wide_gc_mean_weight(str(ref), model) == pytest.approx(0.7)


def test_compute_genome_wide_gc_mean_pure_at_reference(tmp_path):
    """A 100%-AT reference samples every window at GC bin 0 → returns weights[0]."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 2000 + "\n")
    weights = [0.0] * 101
    weights[0] = 0.3
    weights[100] = 1.0
    model = GCBiasModel(weights, window_size=100)
    pysam_index_fasta(ref)
    mean = compute_genome_wide_gc_mean_weight(str(ref), model)
    assert mean == pytest.approx(0.3, abs=1e-6)


def test_compute_genome_wide_gc_mean_mixed_reference(tmp_path):
    """Half-AT/half-GC reference should land between weights[0] and weights[100]."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "A" * 5000 + "G" * 5000 + "\n")
    weights = [0.5] * 101
    weights[100] = 1.0
    model = GCBiasModel(weights, window_size=100)
    pysam_index_fasta(ref)
    mean = compute_genome_wide_gc_mean_weight(str(ref), model)
    # ≈ 50 % windows at bin 0 (weight 0.5) + 50 % at bin 100 (weight 1.0), with a 100 bp
    # transition zone near the midpoint. Should land near 0.75.
    assert 0.7 < mean < 0.8


def pysam_index_fasta(ref_path):
    """Build a .fai for pysam.FastaFile."""
    pysam.faidx(str(ref_path))
