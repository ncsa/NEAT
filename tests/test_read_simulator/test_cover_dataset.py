"""
Tests for sequencing error model in models
"""

import pytest
import numpy as np

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from neat.models import FragmentLengthModel
from neat.variants import SingleNucleotideVariant, Insertion, Deletion
from neat.read_simulator.utils import Options, cover_dataset


def test_cover_dataset():
    """Test that a cover is successfully generated"""
    read_pool = [10] * 2000
    span_length = 100
    target_vector = np.full(100, fill_value=10, dtype=int)
    options = Options(rng_seed=0)
    options.paired_ended = False
    options.read_len = 10
    options.coverage = 10
    fragment_model = FragmentLengthModel(rng=options.rng)

    read1, read2 = cover_dataset(read_pool, span_length, target_vector, options, fragment_model)
    coverage_check = []
    for i in range(span_length):
        # single ended test, only need read1
        cover = [x for x in read1 if i in range(x[0], x[1])]
        coverage_check.append(len(cover))
    assert sum(coverage_check)/len(coverage_check) > 10


def test_paired_cover_dataset():
    """Test that a cover is successfully generated"""
    read_pool = [10] * 2000
    span_length = 100
    target_vector = np.full(100, fill_value=10, dtype=int)
    options = Options(rng_seed=0)
    options.paired_ended = True
    options.read_len = 10
    options.coverage = 10
    fragment_model = FragmentLengthModel(fragment_mean=20, fragment_std=2, fragment_min=10, fragment_max=30, rng=options.rng)

    read1, read2 = cover_dataset(read_pool, span_length, target_vector, options, fragment_model)
    coverage_check = []
    for i in range(span_length):
        # single ended test, only need read1
        cover = [x for x in read1+read2 if i in range(x[0], x[1])]
        coverage_check.append(len(cover))
    assert sum(coverage_check) / len(coverage_check) > 10
