import pytest
import numpy as np

from neat.models import FragmentLengthModel
from neat.read_simulator.utils import Options, cover_dataset


def test_cover_dataset():
    """Test that a cover is successfully generated for different coverage values"""
    read_pool = [10] * 2000
    span_length = 100
    target_vector = np.full(100, fill_value=10, dtype=int)
    options = Options(rng_seed=0)
    options.read_len = 101
    options.paired_ended = True
    options.fragment_mean = 250
    options.fragment_st_dev = 100
    options.output.overwrite_output = True
    fragment_model = FragmentLengthModel(rng=options.rng)

    coverage_values = [1, 2, 5, 10, 25, 50]
    for coverage in coverage_values:
        options.coverage = coverage
        read1, read2 = cover_dataset(read_pool, span_length, target_vector, options, fragment_model)
        coverage_check = []
        for i in range(span_length):
            # single ended test, only need read1
            cover = [x for x in read1 if i in range(x[0], x[1])]
            coverage_check.append(len(cover))
        assert sum(coverage_check)/len(coverage_check) > coverage, f"Coverage check failed for coverage {coverage}"


def test_paired_cover_dataset():
    """Test that a cover is successfully generated for different coverage values"""
    read_pool = [10] * 2000
    span_length = 100
    target_vector = np.full(100, fill_value=10, dtype=int)
    options = Options(rng_seed=0)
    options.read_len = 101
    options.paired_ended = True
    options.fragment_mean = 250
    options.fragment_st_dev = 100
    options.output.overwrite_output = True
    fragment_model = FragmentLengthModel(fragment_mean=20, fragment_std=2, fragment_min=10, fragment_max=30, rng=options.rng)

    coverage_values = [1, 2, 5, 10, 25, 50]
    for coverage in coverage_values:
        options.coverage = coverage
        read1, read2 = cover_dataset(read_pool, span_length, target_vector, options, fragment_model)
        coverage_check = []
        for i in range(span_length):
            # paired ended test, need both read1 and read2
            cover = [x for x in read1 + read2 if i in range(x[0], x[1])]
            coverage_check.append(len(cover))
        assert sum(coverage_check) / len(coverage_check) > coverage, f"Coverage check failed for coverage {coverage}"


def test_various_read_lengths():
    """Test cover_dataset with various read lengths to ensure no errors"""
    read_pool = [10] * 2000
    span_length = 100
    target_vector = np.full(100, fill_value=10, dtype=int)
    options = Options(rng_seed=0)
    options.paired_ended = True
    options.coverage = 10
    options.fragment_mean = 250
    options.fragment_st_dev = 100
    options.output.overwrite_output = True
    fragment_model = FragmentLengthModel(rng=options.rng)

    for read_len in range(10, 251, 10):
        options.read_len = read_len
        try:
            read1, _ = cover_dataset(read_pool, span_length, target_vector, options, fragment_model)
        except Exception as e:
            pytest.fail(f"Test failed for read_len={read_len} with exception: {e}")


def test_fragment_mean_st_dev_combinations():
    """Test cover_dataset with combinations of fragment mean and standard deviation to ensure no errors"""
    read_pool = [10] * 2000
    span_length = 100
    target_vector = np.full(100, fill_value=10, dtype=int)
    options = Options(rng_seed=0)
    options.paired_ended = True
    options.read_len = 50
    options.coverage = 10
    options.output.overwrite_output = True

    fragment_means = [1, 2, 5, 10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000]
    fragment_st_devs = [1, 2, 5, 10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000]

    for mean in fragment_means:
        for st_dev in fragment_st_devs:
            options.fragment_mean = mean
            options.fragment_st_dev = st_dev
            fragment_model = FragmentLengthModel(fragment_mean=mean, fragment_std=st_dev, rng=options.rng)
            try:
                read1, _ = cover_dataset(read_pool, span_length, target_vector, options, fragment_model)
            except Exception as e:
                pytest.fail(f"Test failed for mean={mean}, st_dev={st_dev} with exception: {e}")
