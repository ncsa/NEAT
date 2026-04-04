import pytest

from neat.models import FragmentLengthModel
from neat.read_simulator.utils import Options
from neat.read_simulator.utils.generate_reads import *
from neat.read_simulator.utils.read import *


def _span(a, b):
    """Closed-open span for coverage checks, independent of strand"""
    a = int(a); b = int(b)
    lo, hi = (a, b) if a <= b else (b, a)
    return range(lo, hi)


def _expected_avg_cov_single(requested, span, read_len, frag_mean):
    # Fragment-driven coverage, with finite-span correction
    return requested * (read_len / frag_mean) * max(0.0, (span - int(frag_mean) + 1) / span)


def _expected_avg_cov(requested, span, read_len):
    # Ideal expected average coverage per base for single-end on a finite span, when starts are uniform in [0, span - read_len].
    return requested * max(0, (span - read_len + 1) / span)


def _compute_avg_coverage(reads, span_length, paired):
    """Compute average per-position coverage across the span."""
    cov = [0] * span_length
    for read in reads:
        for pos in range(max(0, read[0]), min(span_length, read[1])):
            cov[pos] += 1
        if paired:
            for pos in range(max(0, read[2]), min(span_length, read[3])):
                cov[pos] += 1
    return sum(cov) / span_length


def test_cover_dataset():
    """Test that a cover is successfully generated for different coverage values"""
    span_length = 5000
    options = Options(rng_seed=0)
    options.read_len = 100
    options.paired_ended = False
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(150, 30)

    coverage_values = [1, 2, 5]
    for coverage in coverage_values:
        options.coverage = coverage
        reads = cover_dataset(span_length, options, fragment_model)
        coverage_check = []
        for i in range(span_length):
            # Single-ended test, only need read1
            cover = [x for x in reads if i in _span(x[0], x[1])]
            coverage_check.append(len(cover))
        print(f"\nRequested Coverage: {coverage}, Average Coverage: {sum(coverage_check) / len(coverage_check)}")
        avg = sum(coverage_check) / len(coverage_check)
        expected = _expected_avg_cov_single(coverage, span_length, options.read_len, fragment_model.fragment_mean)
        assert avg >= 0.85 * expected, f"single-end avg={avg:.3f}, expected≈{expected:.3f} (cov={coverage})"
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
    fragment_model = FragmentLengthModel(300, 30)
    options.fragment_length_model = fragment_model
    options.fragment_mean = fragment_model.fragment_mean
    options.fragment_st_dev = fragment_model.fragment_st_dev

    coverage_values = [1, 2, 5]
    prev_n_reads = None

    for coverage in coverage_values:
        options.coverage = coverage
        reads = cover_dataset(span_length, options, fragment_model)

        expected_pairs = coverage * (span_length / fragment_model.fragment_mean)
        expected_reads = 2.0 * expected_pairs
        n_reads = len(reads)
        assert n_reads >= 0.6 * expected_reads, (
            f"paired-end n_reads={n_reads}, expected≈{expected_reads:.1f} (cov={coverage})"
        )

        assert isinstance(reads, list)
        for read in reads:
            assert len(read) == 4
            assert read[1] >= read[0]
            assert read[3] >= read[2]
            if read[1] - read[0] < 10 or read[3] - read[2] < 10:
                raise AssertionError("failed to filter out a small read length")

        assert len(reads) >= (100 * coverage)/20

        if prev_n_reads is not None:
            assert n_reads >= prev_n_reads
        prev_n_reads = n_reads


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
        reads = cover_dataset(span_length, options, fragment_model)
        assert isinstance(reads, list)


def test_coverage_ploidy_combinations():
    """Exercise cover_dataset over coverage and ploidy combinations to ensure no errors"""
    span_length = 20000
    options = Options(rng_seed=0)
    options.paired_ended = True
    options.read_len = 100
    options.overwrite_output = True
    fragment_model = FragmentLengthModel(300, 100)
    options.fragment_length_model = fragment_model
    options.fragment_mean = fragment_model.fragment_mean
    options.fragment_st_dev = fragment_model.fragment_st_dev

    coverage_values = [1, 2, 5, 10]
    ploidy_values = [1, 2, 4]

    for ploidy in ploidy_values:
        options.ploidy = ploidy
        prev_n_reads = None
        for coverage in coverage_values:
            options.coverage = coverage
            reads = cover_dataset(span_length, options, fragment_model)

            assert isinstance(reads, list)
            for read in reads:
                assert len(read) == 4

            n_reads = len(reads)
            if prev_n_reads is not None:
                assert n_reads >= prev_n_reads
            prev_n_reads = n_reads


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
    fragment_model = FragmentLengthModel(40, 10)

    reads = cover_dataset(span_length, options, fragment_model)
    coverage_check = []
    for i in range(span_length):
        # Single-ended test, only need read1
        cover = [x for x in reads if i in range(x[0], x[1])]
        coverage_check.append(len(cover))
    avg = sum(coverage_check) / len(coverage_check)
    expected = _expected_avg_cov(options.coverage, span_length, options.read_len)
    assert avg >= 0.9 * expected, f"got {avg:.3f}, expected ~{expected:.3f}"


@pytest.mark.parametrize("fragment_mean,target_coverage", [
    (200, 5),
    (200, 10),
    (200, 20),
    (300, 5),
    (300, 10),
    (300, 20),
    (500, 10),
    (500, 20),
])
def test_single_ended_coverage_accuracy(fragment_mean, target_coverage):
    """Single-ended coverage should be within 10% of the requested target."""
    span_length = 10_000
    read_len = 100
    fragment_model = FragmentLengthModel(fragment_mean, 30)

    options = Options(rng_seed=42)
    options.read_len = read_len
    options.paired_ended = False
    options.coverage = target_coverage
    options.overwrite_output = True

    reads = cover_dataset(span_length, options, fragment_model)
    avg = _compute_avg_coverage(reads, span_length, paired=False)

    assert abs(avg - target_coverage) / target_coverage < 0.10, (
        f"Single-ended (frag_mean={fragment_mean}) average coverage {avg:.2f}x "
        f"is more than 10% off target {target_coverage}x"
    )


@pytest.mark.parametrize("fragment_mean,target_coverage", [
    (300, 5),
    (300, 10),
    (300, 20),
    (500, 5),
    (500, 10),
    (500, 20),
    (800, 10),
    (800, 20),
])
def test_paired_ended_coverage_accuracy(fragment_mean, target_coverage):
    """Paired-ended coverage should be within 10% of the requested target."""
    span_length = 10_000
    read_len = 100
    fragment_model = FragmentLengthModel(fragment_mean, 30)

    options = Options(rng_seed=42)
    options.read_len = read_len
    options.paired_ended = True
    options.coverage = target_coverage
    options.overwrite_output = True

    reads = cover_dataset(span_length, options, fragment_model)
    avg = _compute_avg_coverage(reads, span_length, paired=True)

    assert abs(avg - target_coverage) / target_coverage < 0.10, (
        f"Paired-ended (frag_mean={fragment_mean}) average coverage {avg:.2f}x "
        f"is more than 10% off target {target_coverage}x"
    )


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


def test_cigar():
    cigar = ["M"] * 150
    reference = "A" * 150
    read = Read(
        "test",
        (0, 150, 150, 300),
        Seq(reference),
        "chrom",
        0,
        0,
        150,
        15,
        150,
        False,
        True
    )

    cig_str = read.tally_cigar_list(cigar)
    assert cig_str == "150M"

    cigar[11] = "I"
    cig_str = read.tally_cigar_list(cigar)
    assert cig_str == "11M1I138M"

    cigar[137] = "D"
    cig_str = read.tally_cigar_list(cigar)
    assert cig_str == "11M1I125M1D12M"