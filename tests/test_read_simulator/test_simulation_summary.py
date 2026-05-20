"""
Tests for neat/read_simulator/utils/simulation_summary.py

Covers the write_simulation_summary() helper that emits the per-run JSON manifest
consumed by `neat compare-vcfs`, plus the private count helpers.
"""
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from neat.read_simulator.utils.simulation_summary import (
    SCHEMA_VERSION,
    _abs_or_none,
    _count_reads,
    _count_variants,
    _iso_utc,
    write_simulation_summary,
)


# ---------------------------------------------------------------------------
# Synthetic-file builders
# ---------------------------------------------------------------------------

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##contig=<ID=chr2,length=1000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _make_vcf(tmp_path: Path, records: list[tuple[str, int]]) -> Path:
    """records is a list of (chrom, pos) tuples; alt is always 'T' over ref 'A'."""
    path = tmp_path / "variants.vcf"
    with open(path, "w") as fh:
        fh.write(_VCF_HEADER)
        for chrom, pos in records:
            fh.write(f"{chrom}\t{pos}\t.\tA\tT\t.\t.\t.\n")
    return path


def _make_bam(tmp_path: Path, n_reads: int) -> Path:
    """Write a tiny BAM with n_reads identical aligned reads on chr1."""
    path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bf:
        for i in range(n_reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "ACGT"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 100 + i
            a.mapping_quality = 60
            a.cigar = ((0, 4),)
            a.query_qualities = pysam.qualitystring_to_array("IIII")
            bf.write(a)
    return path


def _make_fastq(tmp_path: Path, n_records: int, name: str = "r1.fastq") -> Path:
    path = tmp_path / name
    with open(path, "w") as fh:
        for i in range(n_records):
            fh.write(f"@read{i}\nACGT\n+\nIIII\n")
    return path


def _make_options(**overrides) -> SimpleNamespace:
    """Build a stub Options that exposes the attrs write_simulation_summary reads."""
    defaults = dict(
        reference=None,
        coverage=30,
        read_len=150,
        paired_ended=False,
        fragment_mean=None, fragment_st_dev=None,
        ploidy=2,
        rng_seed=42,
        threads=4,
        mutation_rate=None,
        mutation_bed=None, target_bed=None, discard_bed=None, include_vcf=None,
        mutation_model=None, gc_model=None, error_model=None, fragment_model=None,
        fq1=None, fq2=None, bam=None, vcf=None,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


# ===========================================================================
# _abs_or_none
# ===========================================================================

def test_abs_or_none_returns_none_for_none():
    assert _abs_or_none(None) is None


def test_abs_or_none_resolves_relative_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = _abs_or_none("foo/bar.txt")
    assert Path(result).is_absolute()
    assert result.endswith("foo/bar.txt")


def test_abs_or_none_preserves_absolute_path(tmp_path):
    target = tmp_path / "x.txt"
    assert _abs_or_none(target) == str(target.resolve())


# ===========================================================================
# _iso_utc
# ===========================================================================

def test_iso_utc_format_ends_in_Z():
    """Schema requires an ISO-8601 UTC timestamp with the 'Z' suffix."""
    s = _iso_utc(0)
    assert s.endswith("Z")
    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$", s)


def test_iso_utc_epoch_zero_is_1970():
    assert _iso_utc(0) == "1970-01-01T00:00:00Z"


def test_iso_utc_round_trips_through_datetime():
    epoch = 1_700_000_000.0
    parsed = datetime.fromisoformat(_iso_utc(epoch).replace("Z", "+00:00"))
    assert parsed == datetime.fromtimestamp(epoch, tz=timezone.utc).replace(microsecond=0)


# ===========================================================================
# _count_variants
# ===========================================================================

def test_count_variants_returns_none_pair_for_none_path():
    assert _count_variants(None) == (None, None)


def test_count_variants_returns_none_pair_for_missing_file(tmp_path):
    assert _count_variants(str(tmp_path / "does_not_exist.vcf")) == (None, None)


def test_count_variants_counts_and_groups_by_contig(tmp_path):
    vcf = _make_vcf(tmp_path, [("chr1", 100), ("chr1", 200), ("chr2", 50)])
    total, by_contig = _count_variants(str(vcf))
    assert total == 3
    assert by_contig == {"chr1": 2, "chr2": 1}


def test_count_variants_empty_vcf_returns_zero(tmp_path):
    vcf = _make_vcf(tmp_path, [])
    total, by_contig = _count_variants(str(vcf))
    assert total == 0
    assert by_contig == {}


def test_count_variants_malformed_vcf_returns_none_pair(tmp_path, caplog):
    bad = tmp_path / "bad.vcf"
    bad.write_text("not a real vcf at all\n")
    result = _count_variants(str(bad))
    assert result == (None, None)


# ===========================================================================
# _count_reads
# ===========================================================================

def test_count_reads_returns_none_when_no_inputs():
    assert _count_reads(None, None, None, paired_ended=False) is None


def test_count_reads_prefers_bam_over_fastq(tmp_path):
    bam = _make_bam(tmp_path, 7)
    fq = _make_fastq(tmp_path, 999)
    assert _count_reads(str(bam), str(fq), None, paired_ended=False) == 7


def test_count_reads_falls_back_to_fastq_when_no_bam(tmp_path):
    fq = _make_fastq(tmp_path, 12)
    assert _count_reads(None, str(fq), None, paired_ended=False) == 12


def test_count_reads_doubles_for_paired_ended_fastq(tmp_path):
    fq = _make_fastq(tmp_path, 10, name="r1.fastq")
    assert _count_reads(None, str(fq), None, paired_ended=True) == 20


def test_count_reads_returns_none_for_missing_bam_and_missing_fq(tmp_path):
    assert _count_reads(
        str(tmp_path / "nope.bam"), str(tmp_path / "nope.fq"), None, paired_ended=False
    ) is None


# ===========================================================================
# write_simulation_summary — schema shape and content
# ===========================================================================

def test_write_simulation_summary_writes_file_in_output_dir(tmp_path):
    options = _make_options()
    path = write_simulation_summary(
        options=options, output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0,
        contigs_simulated=[],
    )
    assert path == tmp_path / "simulation_summary.json"
    assert path.is_file()


def test_write_simulation_summary_has_expected_top_level_keys(tmp_path):
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    summary = json.loads((tmp_path / "simulation_summary.json").read_text())
    assert set(summary) == {"schema_version", "neat_version", "run", "config", "outputs", "delivered"}


def test_write_simulation_summary_schema_version_matches_constant(tmp_path):
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    summary = json.loads((tmp_path / "simulation_summary.json").read_text())
    assert summary["schema_version"] == SCHEMA_VERSION


def test_write_simulation_summary_neat_version_is_populated(tmp_path):
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    summary = json.loads((tmp_path / "simulation_summary.json").read_text())
    assert isinstance(summary["neat_version"], str) and summary["neat_version"]


def test_write_simulation_summary_echoes_config_fields(tmp_path):
    options = _make_options(coverage=42, read_len=200, paired_ended=True, ploidy=4, rng_seed=99, threads=8)
    write_simulation_summary(
        options=options, output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    cfg = json.loads((tmp_path / "simulation_summary.json").read_text())["config"]
    assert cfg["coverage"] == 42
    assert cfg["read_len"] == 200
    assert cfg["paired_ended"] is True
    assert cfg["ploidy"] == 4
    assert cfg["rng_seed"] == 99
    assert cfg["threads"] == 8


def test_write_simulation_summary_optional_config_fields_default_to_null(tmp_path):
    """All Path-typed optional fields and rate must serialize as JSON null when absent."""
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    cfg = json.loads((tmp_path / "simulation_summary.json").read_text())["config"]
    optional_keys = [
        "reference", "mutation_rate", "mutation_bed", "target_bed", "discard_bed",
        "include_vcf", "mutation_model", "gc_model", "error_model", "fragment_model",
    ]
    for k in optional_keys:
        assert cfg[k] is None, f"expected {k} to be null"


def test_write_simulation_summary_resolves_paths_to_absolute(tmp_path, monkeypatch):
    """Relative inputs (cwd, reference, config_path) must serialize as absolute paths."""
    monkeypatch.chdir(tmp_path)
    ref = tmp_path / "ref.fa"
    ref.touch()
    out_rel = tmp_path / "out_rel"
    out_rel.mkdir()  # runner creates output_dir before calling the helper
    options = _make_options(reference="ref.fa")  # relative
    write_simulation_summary(
        options=options, output_dir="out_rel",
        file_prefix="run", config_path="cfg.yml", analysis_start=0.0,
        contigs_simulated=[],
    )
    summary = json.loads((out_rel / "simulation_summary.json").read_text())
    assert Path(summary["run"]["output_dir"]).is_absolute()
    assert Path(summary["run"]["config_file"]).is_absolute()
    assert Path(summary["config"]["reference"]).is_absolute()


def test_write_simulation_summary_contigs_simulated_preserves_order(tmp_path):
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0,
        contigs_simulated=["chr3", "chr1", "chr2"],
    )
    summary = json.loads((tmp_path / "simulation_summary.json").read_text())
    assert summary["delivered"]["contigs_simulated"] == ["chr3", "chr1", "chr2"]


def test_write_simulation_summary_records_reference_contigs(tmp_path):
    """`reference_contigs` captures the full FASTA contig set, separate from
    `contigs_simulated` (which is what the simulator iterated over)."""
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0,
        contigs_simulated=["chr1"],
        reference_contigs=["chr1", "chr2", "chrX", "chrM"],
    )
    delivered = json.loads((tmp_path / "simulation_summary.json").read_text())["delivered"]
    assert delivered["reference_contigs"] == ["chr1", "chr2", "chrX", "chrM"]
    assert delivered["contigs_simulated"] == ["chr1"]


def test_write_simulation_summary_reference_contigs_defaults_to_contigs_simulated(tmp_path):
    """When reference_contigs is omitted, fall back to contigs_simulated so old
    callers that don't pass it still get a populated field."""
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0,
        contigs_simulated=["chr1", "chr2"],
    )
    delivered = json.loads((tmp_path / "simulation_summary.json").read_text())["delivered"]
    assert delivered["reference_contigs"] == ["chr1", "chr2"]


def test_write_simulation_summary_outputs_none_when_nothing_produced(tmp_path):
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    outputs = json.loads((tmp_path / "simulation_summary.json").read_text())["outputs"]
    assert outputs == {"fastq": None, "bam": None, "vcf": None}


def test_write_simulation_summary_outputs_populated_for_multi_output_run(tmp_path):
    vcf = _make_vcf(tmp_path, [("chr1", 100), ("chr1", 200), ("chr2", 50)])
    bam = _make_bam(tmp_path, 5)
    fq1 = _make_fastq(tmp_path, 3, "r1.fastq")
    fq2 = _make_fastq(tmp_path, 3, "r2.fastq")
    options = _make_options(paired_ended=True, vcf=vcf, bam=bam, fq1=fq1, fq2=fq2)
    write_simulation_summary(
        options=options, output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0,
        contigs_simulated=["chr1", "chr2"],
    )
    summary = json.loads((tmp_path / "simulation_summary.json").read_text())
    assert summary["outputs"]["vcf"] == str(vcf.resolve())
    assert summary["outputs"]["bam"] == str(bam.resolve())
    assert summary["outputs"]["fastq"] == [str(fq1.resolve()), str(fq2.resolve())]
    assert summary["delivered"]["total_variants"] == 3
    assert summary["delivered"]["variants_by_contig"] == {"chr1": 2, "chr2": 1}
    assert summary["delivered"]["total_reads"] == 5  # BAM preferred over FASTQ


def test_write_simulation_summary_duration_matches_elapsed(tmp_path):
    """duration_seconds should be the difference between now and analysis_start."""
    import time
    start = time.time() - 10  # ten seconds ago
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=start, contigs_simulated=[],
    )
    summary = json.loads((tmp_path / "simulation_summary.json").read_text())
    assert 9.5 <= summary["run"]["duration_seconds"] <= 11.0


def test_write_simulation_summary_no_tmp_file_remains(tmp_path):
    """Atomic write contract: .json.tmp must not linger after success."""
    write_simulation_summary(
        options=_make_options(), output_dir=tmp_path, file_prefix="run",
        config_path=tmp_path / "cfg.yml", analysis_start=0.0, contigs_simulated=[],
    )
    assert not (tmp_path / "simulation_summary.json.tmp").exists()
    assert (tmp_path / "simulation_summary.json").exists()
