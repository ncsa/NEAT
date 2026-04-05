"""
Tests for neat/read_simulator/runner.py

Unit tests cover filter_thread_variants and filter_bed_regions.
Integration test exercises read_simulator_runner end-to-end with a tiny
reference and a minimal config, producing only FASTQ output so that
pysam/bcftools post-processing is not triggered.
"""
import gzip
import textwrap
from pathlib import Path

import numpy as np
import pytest

from neat.read_simulator.runner import (
    filter_thread_variants,
    filter_bed_regions,
    read_simulator_runner,
)
from neat.variants import ContigVariants, SingleNucleotideVariant, Deletion, Insertion


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def contig_variants():
    """ContigVariants with SNVs at positions 10, 50, 100, 200."""
    cv = ContigVariants()
    for pos in (10, 50, 100, 200):
        cv.add_variant(SingleNucleotideVariant(pos, "T", np.array([0, 1]), "42"))
    return cv


# ===========================================================================
# filter_thread_variants
# ===========================================================================

def test_filter_thread_variants_all_in_range(contig_variants):
    result = filter_thread_variants(contig_variants, (0, 300))
    assert sorted(result.variant_locations) == [10, 50, 100, 200]


def test_filter_thread_variants_none_in_range(contig_variants):
    result = filter_thread_variants(contig_variants, (500, 1000))
    assert result.variant_locations == []


def test_filter_thread_variants_partial_range(contig_variants):
    result = filter_thread_variants(contig_variants, (50, 150))
    assert sorted(result.variant_locations) == [50, 100]


def test_filter_thread_variants_lower_bound_inclusive(contig_variants):
    """coords[0] == variant position → included."""
    result = filter_thread_variants(contig_variants, (10, 50))
    assert 10 in result.variant_locations


def test_filter_thread_variants_upper_bound_exclusive(contig_variants):
    """coords[1] == variant position → excluded."""
    result = filter_thread_variants(contig_variants, (10, 50))
    assert 50 not in result.variant_locations


def test_filter_thread_variants_empty_input():
    result = filter_thread_variants(ContigVariants(), (0, 1000))
    assert result.variant_locations == []


def test_filter_thread_variants_returns_contig_variants(contig_variants):
    result = filter_thread_variants(contig_variants, (0, 300))
    assert isinstance(result, ContigVariants)


def test_filter_thread_variants_preserves_variant_data(contig_variants):
    """Filtered variants should retain their original objects."""
    result = filter_thread_variants(contig_variants, (0, 60))
    for loc in result.variant_locations:
        variants = result.contig_variants[loc]
        assert len(variants) == 1
        assert isinstance(variants[0], SingleNucleotideVariant)


# ===========================================================================
# filter_bed_regions
# ===========================================================================

_REGIONS = [
    (0,   200,  0.01),
    (200, 400,  0.02),
    (400, 600,  0.03),
    (600, 1000, 0.04),
]


def test_filter_bed_regions_block_within_one_region():
    result = filter_bed_regions(_REGIONS, (50, 150))
    assert result == [(0, 200, 0.01)]


def test_filter_bed_regions_block_spans_two_regions():
    result = filter_bed_regions(_REGIONS, (150, 250))
    assert (0, 200, 0.01) in result
    assert (200, 400, 0.02) in result


def test_filter_bed_regions_block_spans_all():
    result = filter_bed_regions(_REGIONS, (0, 1000))
    assert len(result) == len(_REGIONS)


def test_filter_bed_regions_block_outside_all():
    result = filter_bed_regions(_REGIONS, (1100, 1500))
    assert result == []


def test_filter_bed_regions_empty_regions():
    result = filter_bed_regions([], (0, 1000))
    assert result == []


def test_filter_bed_regions_block_touching_region_start():
    """Block ends exactly at region boundary → the touching region is included."""
    # coords (200, 400) → region (200, 400) shares its left edge with coords[0]
    result = filter_bed_regions(_REGIONS, (200, 400))
    assert (200, 400, 0.02) in result


def test_filter_bed_regions_single_region_fully_inside_block():
    """Region fully contained within block coordinates."""
    result = filter_bed_regions([(300, 350, 0.05)], (200, 400))
    assert result == [(300, 350, 0.05)]


def test_filter_bed_regions_returns_list():
    result = filter_bed_regions(_REGIONS, (50, 150))
    assert isinstance(result, list)


# ===========================================================================
# Integration test — read_simulator_runner (FASTQ output only)
# ===========================================================================

def _write_ref(path: Path, seq: str = "ACGT" * 100) -> Path:
    """Write a minimal single-contig FASTA reference."""
    path.write_text(f">chr1\n{seq}\n", encoding="utf-8")
    return path


def _write_config(path: Path, ref_path: Path, **overrides) -> Path:
    defaults = {
        "reference": str(ref_path),
        "produce_fastq": "true",
        "produce_vcf": "false",
        "produce_bam": "false",
        "read_len": 50,
        "coverage": 2,
        "rng_seed": 42,
        "overwrite_output": "true",
        "cleanup_splits": "true",
    }
    defaults.update(overrides)
    lines = "\n".join(f"{k}: {v}" for k, v in defaults.items())
    path.write_text(lines + "\n", encoding="utf-8")
    return path


def test_runner_produces_fastq_output(tmp_path):
    """
    End-to-end: runner with a 400 bp reference, low coverage, FASTQ only.
    Verifies fq1 is created and contains at least one FASTQ record.
    """
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(tmp_path / "conf.yml", ref)
    out_dir = tmp_path / "out"

    read_simulator_runner(str(cfg), str(out_dir), "test")

    fq_files = list(out_dir.glob("*.fastq.gz"))
    assert len(fq_files) >= 1, "Expected at least one FASTQ output file"

    # Verify the file contains valid FASTQ content
    with gzip.open(fq_files[0], "rt") as fh:
        first_line = fh.readline()
    assert first_line.startswith("@"), "FASTQ file should start with '@'"


def test_runner_paired_end_produces_two_fastqs(tmp_path):
    """Paired-end mode produces both fq1 and fq2."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        paired_ended="true",
        fragment_mean=150,
        fragment_st_dev=25,
    )
    out_dir = tmp_path / "out"

    read_simulator_runner(str(cfg), str(out_dir), "test")

    fq_files = sorted(out_dir.glob("*.fastq.gz"))
    assert len(fq_files) == 2, f"Expected 2 FASTQ files, found {len(fq_files)}"


def test_runner_creates_output_dir_if_missing(tmp_path):
    """Output directory is created if it doesn't exist."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(tmp_path / "conf.yml", ref)
    out_dir = tmp_path / "nested" / "output"

    assert not out_dir.exists()
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert out_dir.is_dir()


def test_runner_output_prefix_applied(tmp_path):
    """Output files use the supplied prefix."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(tmp_path / "conf.yml", ref)
    out_dir = tmp_path / "out"

    read_simulator_runner(str(cfg), str(out_dir), "myprefix")

    output_files = list(out_dir.glob("myprefix*"))
    assert len(output_files) >= 1, "No files found with the expected prefix"


def test_runner_reproducible_with_same_seed(tmp_path):
    """Two runs with the same seed produce identical FASTQ output."""
    ref = _write_ref(tmp_path / "ref.fa")

    out1 = tmp_path / "run1"
    cfg1 = _write_config(tmp_path / "conf1.yml", ref, rng_seed=7)
    read_simulator_runner(str(cfg1), str(out1), "rep")

    out2 = tmp_path / "run2"
    cfg2 = _write_config(tmp_path / "conf2.yml", ref, rng_seed=7)
    read_simulator_runner(str(cfg2), str(out2), "rep")

    fq1_a = sorted(out1.glob("*.fastq.gz"))[0]
    fq1_b = sorted(out2.glob("*.fastq.gz"))[0]

    with gzip.open(fq1_a, "rt") as a, gzip.open(fq1_b, "rt") as b:
        assert a.read() == b.read(), "Seeded runs should be identical"


def test_runner_with_vcf_output(tmp_path):
    """Runner with produce_vcf=true creates a VCF output file."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        produce_vcf="true",
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"

    read_simulator_runner(str(cfg), str(out_dir), "test")

    vcf_files = list(out_dir.glob("*.vcf.gz"))
    assert len(vcf_files) == 1, "Expected exactly one VCF output file"
    assert vcf_files[0].stat().st_size > 0
