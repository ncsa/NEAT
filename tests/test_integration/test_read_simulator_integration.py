from pathlib import Path

import gzip
import pytest

from neat.read_simulator.runner import read_simulator_runner


DATA_DIR = Path(__file__).resolve().parents[2] / "data"


@pytest.mark.integration
def test_single_end_quick_run_h1n1(tmp_path: Path):
    # Use a minimal config derived from H1N1 that keeps runtime small
    cfg = tmp_path / "cfg_se.yml"
    cfg.write_text(
        "\n".join([
            "reference: data/H1N1.fa",
            "read_len: 50",
            "coverage: 1",
            "paired_ended: false",
            "produce_fastq: true",
            "produce_bam: false",
            "produce_vcf: false",
            "overwrite_output: true",
            "threads: 1",
        ]),
        encoding="utf-8",
    )

    outdir = tmp_path / "out"
    prefix = "itest_se"

    read_simulator_runner(str(cfg), str(outdir), prefix)

    # Should produce a single fastq.gz
    fq = outdir / f"{prefix}.fastq.gz"
    assert fq.exists() and fq.stat().st_size > 0

    # Optionally inspect first few lines to ensure text
    with gzip.open(fq, "rt", encoding="utf-8") as fh:
        head = [next(fh) for _ in range(4)]
    assert head[0].startswith("@")


@pytest.mark.integration
def test_paired_end_quick_run_h1n1(tmp_path: Path):
    # Minimal paired-end config; small coverage to keep test fast
    cfg = tmp_path / "cfg_pe.yml"
    cfg.write_text(
        "\n".join([
            "reference: data/H1N1.fa",
            "read_len: 75",
            "coverage: 1",
            "paired_ended: true",
            "fragment_mean: 250",
            "fragment_st_dev: 50",
            "produce_fastq: true",
            "produce_bam: false",
            "produce_vcf: false",
            "overwrite_output: true",
            "threads: 1",
        ]),
        encoding="utf-8",
    )

    outdir = tmp_path / "out"
    prefix = "itest_pe"

    read_simulator_runner(str(cfg), str(outdir), prefix)

    fq1 = outdir / f"{prefix}_r1.fastq.gz"
    fq2 = outdir / f"{prefix}_r2.fastq.gz"
    assert fq1.exists() and fq1.stat().st_size > 0
    assert fq2.exists() and fq2.stat().st_size > 0

    # Spot-check content structure
    with gzip.open(fq1, "rt", encoding="utf-8") as fh:
        l1 = next(fh)
    assert l1.startswith("@")
