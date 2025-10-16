from math import ceil

import pytest
import numpy as np

from neat.models import FragmentLengthModel
from neat.read_simulator.utils import Options
from neat.read_simulator.utils.generate_reads import *


def test_cover_dataset():
    """Test that a cover is successfully generated for different coverage values"""
    span_length = 100
    options = Options(rng_seed=0)
    options.read_len = 10
    options.paired_ended = False
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(15, 5)

    coverage_values = [1, 2, 5]
    for coverage in coverage_values:
        options.coverage = coverage
        reads = cover_dataset(span_length, options, fragment_model)
        coverage_check = []
        for i in range(span_length):
            # single ended test, only need read1
            cover = [x for x in reads if i in range(x[0], x[1])]
            coverage_check.append(len(cover))
        print(f"\nRequested Coverage: {coverage}, Average Coverage: {sum(coverage_check) / len(coverage_check)}")
        assert sum(coverage_check)/len(coverage_check) > 0.9 * coverage, f"Coverage check failed for coverage {coverage}"
        for read in reads:
            if read[1] - read[0] < 10:
                raise AssertionError("failed to filter out a small read length")
        assert len(reads) >= (100 * coverage)/10


def test_paired_cover_dataset():
    """Test that a cover is successfully generated for different coverage values"""
    span_length = 10000
    options = Options(rng_seed=0)
    options.read_len = 100
    options.paired_ended = True
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(250, 10)

    coverage_values = [1, 2, 5]
    for coverage in coverage_values:
        options.coverage = coverage
        reads = cover_dataset(span_length, options, fragment_model)
        coverage_check = []
        for i in range(span_length):
            # paired ended test, need both read1 and read2
            cover = [x for x in reads if (i in range(x[0], x[1]) or i in range(x[2], x[3]))]
            coverage_check.append(len(cover))
        print(f"\nRequested Coverage: {coverage}, Average Coverage: {sum(coverage_check) / len(coverage_check)}")
        assert sum(coverage_check) / len(coverage_check) > 0.85 * coverage, \
            f"Coverage check failed for coverage {coverage}"
        for read in reads:
            if read[1] - read[0] < 10 or read[3] - read[2] < 10:
                raise AssertionError("failed to filter out a small read length")
        assert len(reads) >= (100 * coverage)/20


def test_various_read_lengths():
    """Test cover_dataset with various read lengths to ensure no errors"""
    span_length = 1000
    options = Options(rng_seed=0)
    options.paired_ended = True
    options.coverage = 10
    options.fragment_mean = 250
    options.fragment_st_dev = 100
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(250, 100)

    for read_len in range(10, 251, 10):
        options.read_len = read_len
        try:
            reads = cover_dataset(span_length, options, fragment_model)
            assert isinstance(reads, list)
        except Exception as e:
            pytest.fail(f"Test failed for read_len={read_len} with exception: {e}")


def test_fragment_mean_st_dev_combinations():
    """Test cover_dataset with combinations of fragment mean and standard deviation to ensure no errors"""
    span_length = 5000
    options = Options(rng_seed=0)
    options.paired_ended = False
    options.read_len = 101
    options.coverage = 2
    options.overwrite_output = True

    fragment_means = [100, 150, 200, 250,]
    fragment_st_devs = [1, 5, 25, 50]

    for mean in fragment_means:
        for st_dev in fragment_st_devs:
            options.fragment_mean = mean
            options.fragment_st_dev = st_dev
            fragment_model = FragmentLengthModel(mean, st_dev)
            frags = fragment_model.generate_fragments(20, options.rng)
            assert len(frags) == 20
            assert fragment_model.fragment_mean == mean
            assert fragment_model.fragment_st_dev == st_dev

def test_coverage_ploidy_combinations():
    """Test cover_dataset with various combinations of coverage and ploidy values to ensure no errors"""
    span_length = 10000
    options = Options(rng_seed=0)
    options.paired_ended = True
    options.read_len = 100
    options.fragment_mean = 250
    options.fragment_st_dev = 100
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(250, 100)

    coverage_values = [1, 2, 5, 10]
    ploidy_values = [1, 2, 4]

    for coverage in coverage_values:
        for ploidy in ploidy_values:
            options.coverage = coverage
            options.ploidy = ploidy
            reads = cover_dataset(span_length, options, fragment_model)
            coverage_check = []
            for i in range(span_length):
                # paired ended test, need either read1 or read2 to cover position
                cover = [x for x in reads if (i in range(x[0], x[1]) or i in range(x[2], x[3]))]
                coverage_check.append(len(cover))
            assert sum(coverage_check) / len(coverage_check) > 0.9 * coverage, (
                f"Coverage check failed for coverage {coverage} and ploidy {ploidy}")

def test_single_ended_mode():
    """Test cover_dataset in single-ended mode for various configurations"""
    span_length = 100
    options = Options(rng_seed=0)
    options.read_len = 50
    options.paired_ended = False
    options.fragment_mean = 250
    options.fragment_st_dev = 100
    options.coverage = 10
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(250, 100)

    try:
        reads = cover_dataset(span_length, options, fragment_model)
        coverage_check = []
        for i in range(span_length):
            # Single-ended test, only need read1
            cover = [x for x in reads if i in range(x[0], x[1])]
            coverage_check.append(len(cover))
        assert sum(coverage_check) / len(coverage_check) > 0.9 * options.coverage, "Coverage check failed in single-ended mode"
    except Exception as e:
        pytest.fail(f"Test failed in single-ended mode with exception: {e}")


def test_overlaps():
    comparsion_interval = (120, 400)

    interval1 = (0, 151)
    assert overlaps(interval1, comparsion_interval)  # overlaps
    interval1 = (385, 421)
    assert overlaps(interval1, comparsion_interval)  # overlaps
    interval1 = (130, 281)
    assert overlaps(interval1, comparsion_interval)  # overlaps
    interval1 = (120, 400)
    assert overlaps(interval1, comparsion_interval)  # overlaps
    interval1 = (100, 500)
    assert overlaps(interval1, comparsion_interval)  # overlaps
    interval1 = (0, 100)
    assert not overlaps(interval1, comparsion_interval)  # doesn't overlap
    interval1 = (0, 120)
    assert not overlaps(interval1, comparsion_interval)  # doesn't overlap
    interval1 = (400, 450)
    assert not overlaps(interval1, comparsion_interval)  # doesn't overlap
    interval1 = (500, 700)
    assert not overlaps(interval1, comparsion_interval)  # doesn't overlap
