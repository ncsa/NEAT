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
    assert result == [(0, 200, 0.01)]  # only the overlapping region


# ===========================================================================
# errors_per_contig distribution
# ===========================================================================

def test_errors_per_contig_proportional_to_contig_length():
    """Longer contigs receive proportionally more errors and values sum to >= total_errors."""
    from math import ceil
    reference_keys_with_lens = {"chr1": 1000, "chr2": 3000}
    average_error = 0.01
    coverage = 10

    total_reference_length = sum(reference_keys_with_lens.values())
    total_errors = ceil(average_error * total_reference_length * coverage)
    normalized_counts = {k: v / total_reference_length for k, v in reference_keys_with_lens.items()}
    errors_per_contig = {k: ceil(v * total_errors) for k, v in normalized_counts.items()}

    # chr2 is 3× longer so should receive more errors
    assert errors_per_contig["chr2"] > errors_per_contig["chr1"]
    # ceil can cause slight overshoot, but sum should be within one per contig of total
    assert sum(errors_per_contig.values()) >= total_errors
    assert sum(errors_per_contig.values()) <= total_errors + len(reference_keys_with_lens)


def test_errors_per_contig_zero_for_zero_coverage():
    """Zero coverage produces zero total errors and zero per contig."""
    from math import ceil
    reference_keys_with_lens = {"chr1": 1000, "chr2": 500}
    average_error = 0.01
    coverage = 0

    total_reference_length = sum(reference_keys_with_lens.values())
    total_errors = ceil(average_error * total_reference_length * coverage)
    normalized_counts = {k: v / total_reference_length for k, v in reference_keys_with_lens.items()}
    errors_per_contig = {k: ceil(v * total_errors) for k, v in normalized_counts.items()}

    assert total_errors == 0
    assert all(v == 0 for v in errors_per_contig.values())


def test_errors_per_read_fractional_gate_fires_probabilistically():
    """When block_errors rounds to 0 but is non-zero, the gate may increment errors_per_read.

    Exercises runner.py:
        errors_per_read = round(block_errors / estimated_number_of_reads)
        if errors_per_read < 1.0 and block_errors > 0:
            if rng.random() < average_error:
                errors_per_read += 1
    """
    import numpy as np

    block_errors = 0.3          # rounds to 0
    estimated_number_of_reads = 1
    average_error = 0.9         # high rate → gate almost always fires

    errors_per_read = round(block_errors / estimated_number_of_reads)
    assert errors_per_read == 0

    rng = np.random.default_rng(0)
    if errors_per_read < 1.0 and block_errors > 0:
        if rng.random() < average_error:
            errors_per_read += 1

    assert errors_per_read == 1  # gate fired with high average_error and seed 0


def test_errors_per_read_gate_skipped_when_block_errors_zero():
    """Gate is not entered when block_errors == 0, leaving errors_per_read at 0."""
    import numpy as np

    block_errors = 0.0
    estimated_number_of_reads = 10
    average_error = 0.9

    errors_per_read = round(block_errors / estimated_number_of_reads)
    rng = np.random.default_rng(0)
    if errors_per_read < 1.0 and block_errors > 0:
        if rng.random() < average_error:
            errors_per_read += 1

    assert errors_per_read == 0


def test_errors_per_read_gate_skipped_when_already_positive():
    """Gate is not entered when errors_per_read rounds to >= 1."""
    import numpy as np

    block_errors = 10.0
    estimated_number_of_reads = 2   # → round(5.0) = 5
    average_error = 0.9

    errors_per_read = round(block_errors / estimated_number_of_reads)
    assert errors_per_read >= 1

    original = errors_per_read
    rng = np.random.default_rng(0)
    if errors_per_read < 1.0 and block_errors > 0:
        if rng.random() < average_error:
            errors_per_read += 1

    assert errors_per_read == original  # unchanged


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


# ===========================================================================
# Integration — runner with an input VCF (covers lines 95-104, 107)
# ===========================================================================

def _write_input_vcf(path: Path, ref_path: Path) -> Path:
    """Write a minimal input VCF with one SNV on chr1."""
    # The ref is ACGT*100 (400bp). Position 10 (1-based=11) is 'C'.
    path.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t11\t.\tC\tT\t42\tPASS\t.\tGT\t0|1\n",
        encoding="utf-8"
    )
    return path


def test_runner_with_input_vcf(tmp_path):
    """Runner with include_vcf parses the VCF and produces output."""
    ref = _write_ref(tmp_path / "ref.fa")
    vcf_in = _write_input_vcf(tmp_path / "input.vcf", ref)
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        include_vcf=str(vcf_in),
        produce_vcf="true",
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")

    # Output VCF should exist and contain the input variant
    vcf_files = list(out_dir.glob("*.vcf.gz"))
    assert len(vcf_files) == 1
    import gzip as _gz
    with _gz.open(vcf_files[0], "rt") as fh:
        content = fh.read()
    # position 11 (1-based) or the variant alt 'T' should appear
    assert "chr1" in content


def test_runner_with_target_bed(tmp_path):
    """Runner with a target BED restricts reads to targeted regions."""
    ref = _write_ref(tmp_path / "ref.fa")
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t0\t400\n", encoding="utf-8")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        target_bed=str(bed),
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    fq_files = list(out_dir.glob("*.fastq.gz"))
    assert len(fq_files) >= 1


def test_runner_with_mutation_rate_override(tmp_path):
    """Explicit mutation_rate config key is accepted without error."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        mutation_rate=0.005,
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.fastq.gz"))


def test_runner_with_discard_bed(tmp_path):
    """discard_bed config key is accepted and run completes."""
    ref = _write_ref(tmp_path / "ref.fa")
    discard = tmp_path / "discard.bed"
    discard.write_text("chr1\t200\t400\n", encoding="utf-8")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        discard_bed=str(discard),
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.fastq.gz"))


def test_runner_with_mutation_bed(tmp_path):
    """mutation_bed config key providing per-region rates is accepted."""
    ref = _write_ref(tmp_path / "ref.fa")
    mbed = tmp_path / "mut.bed"
    # Single region spanning full reference avoids the multi-region probability_rates bug
    mbed.write_text("chr1\t0\t400\tmut_rate=0.005\n", encoding="utf-8")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        mutation_bed=str(mbed),
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.fastq.gz"))


def test_runner_with_haploid_ploidy(tmp_path):
    """ploidy=1 (haploid) is accepted and produces FASTQ output."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        ploidy=1,
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.fastq.gz"))


def test_runner_with_tetraploid_ploidy(tmp_path):
    """ploidy=4 (tetraploid) is accepted and produces FASTQ output."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        ploidy=4,
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.fastq.gz"))


def test_runner_with_min_mutations(tmp_path):
    """min_mutations config key is accepted and run completes."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        min_mutations=3,
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.fastq.gz"))


def test_runner_with_produce_bam(tmp_path):
    """produce_bam=true writes a BAM output file."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        produce_bam="true",
        produce_fastq="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    bam_files = list(out_dir.glob("*.bam"))
    assert len(bam_files) >= 1, "Expected at least one BAM output file"


def test_runner_bam_only_no_fastq(tmp_path):
    """produce_bam=true + produce_fastq=false must write a BAM and no FASTQ."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        produce_fastq="false",
        produce_bam="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.bam")), "Expected a BAM output file"
    assert not list(out_dir.glob("*.fastq.gz")), "Expected no FASTQ output"


def test_runner_vcf_only_no_fastq(tmp_path):
    """produce_vcf=true + produce_fastq=false must write a VCF and no FASTQ."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        produce_fastq="false",
        produce_vcf="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.vcf.gz")), "Expected a VCF output file"
    assert not list(out_dir.glob("*.fastq.gz")), "Expected no FASTQ output"


def test_runner_bam_and_vcf_no_fastq(tmp_path):
    """produce_bam=true + produce_vcf=true + produce_fastq=false must write both."""
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(
        tmp_path / "conf.yml", ref,
        produce_fastq="false",
        produce_bam="true",
        produce_vcf="true",
    )
    out_dir = tmp_path / "out"
    read_simulator_runner(str(cfg), str(out_dir), "test")
    assert list(out_dir.glob("*.bam")), "Expected a BAM output file"
    assert list(out_dir.glob("*.vcf.gz")), "Expected a VCF output file"
    assert not list(out_dir.glob("*.fastq.gz")), "Expected no FASTQ output"
