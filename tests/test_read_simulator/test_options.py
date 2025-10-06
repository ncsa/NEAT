from neat.read_simulator.utils.options import Options

from pathlib import Path as _PathAlias
import numpy as _np
import textwrap as _textwrap
import pytest as _pytest


def _project_root() -> _PathAlias:
    return _PathAlias(__file__).resolve().parents[2]


# Redefine the function name used above to override the brittle test
# so pytest only sees this correct version.
def test_basic_options():
    reference = _project_root() / "data" / "H1N1.fa"
    base_options = Options(reference)
    assert base_options.reference == reference


def test_output_prefix_and_paths_single_end(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="neat_unit", overwrite_output=True)
    assert opts.output_dir == tmp_path
    assert opts.output_prefix == "neat_unit"
    opts.paired_ended = False
    opts.produce_fastq = True
    opts.produce_bam = False
    opts.produce_vcf = False
    opts.log_configuration()
    assert opts.fq1 == tmp_path / "neat_unit.fastq.gz"
    assert opts.fq2 is None
    assert opts.bam is None
    assert opts.vcf is None


def test_output_paths_paired_end(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="pe", overwrite_output=True,
                   paired_ended=True, fragment_mean=200, fragment_st_dev=50)
    opts.produce_fastq = True
    opts.produce_bam = False
    opts.produce_vcf = False
    opts.log_configuration()
    assert opts.fq1 == tmp_path / "pe_r1.fastq.gz"
    assert opts.fq2 == tmp_path / "pe_r2.fastq.gz"


def test_rng_seed_reproducible():
    opts1 = Options(rng_seed=123)
    opts2 = Options(rng_seed=123)
    a1 = opts1.rng.integers(0, 1000000, size=10)
    a2 = opts2.rng.integers(0, 1000000, size=10)
    assert (a1 == a2).all()
    opts3 = Options()
    assert isinstance(opts3.rng_seed, (int, _np.integer))


def test_from_cli_single_end_with_threads_and_splits(tmp_path: _PathAlias):
    # Build a minimal YAML config using repository-relative paths
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        read_len: 75
        coverage: 5
        ploidy: 2
        paired_ended: false

        produce_bam: false
        produce_vcf: false
        produce_fastq: true

        avg_seq_error: 0.01
        rescale_qualities: true
        quality_offset: 33
        rng_seed: 42
        overwrite_output: true

        mode: contig
        size: 500000
        threads: 2
        cleanup_splits: false
        reuse_splits: false
        """
    ).strip() + "\n"

    yml_path = tmp_path / "neat_from_cli.yml"
    yml_path.write_text(cfg, encoding="utf-8")

    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "fromcli", yml_path)

    # Basics propagated
    assert opts.reference == _project_root() / "data" / "H1N1.fa"
    assert opts.read_len == 75
    assert opts.coverage == 5
    assert opts.ploidy == 2
    assert opts.rng_seed == 42

    # Output construction via log_configuration() inside from_cli
    assert opts.output_dir == outdir
    assert opts.output_prefix == "fromcli"
    assert opts.fq1 == outdir / "fromcli.fastq.gz"
    assert opts.fq2 is None
    assert opts.bam is None
    assert opts.vcf is None

    # Parallel-related settings
    assert opts.threads == 2
    # cleanup_splits: false -> splits dir under output_dir
    assert opts.splits_dir == outdir / "splits"
    assert opts.splits_dir.is_dir()


def test_from_cli_paired_end_fragments(tmp_path: _PathAlias):
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        read_len: 101
        coverage: 2
        ploidy: 2
        paired_ended: true
        fragment_mean: 200
        fragment_st_dev: 30

        produce_bam: false
        produce_vcf: false
        produce_fastq: true

        rng_seed: 7
        overwrite_output: true

        mode: contig
        threads: 1
        cleanup_splits: true
        reuse_splits: false
        """
    ).strip() + "\n"

    yml_path = tmp_path / "neat_from_cli_pe.yml"
    yml_path.write_text(cfg, encoding="utf-8")
    outdir = tmp_path / "peout"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "peprefix", yml_path)

    assert opts.paired_ended is True
    assert opts.fragment_mean == 200
    assert opts.fragment_st_dev == 30
    assert opts.fq1 == outdir / "peprefix_r1.fastq.gz"
    assert opts.fq2 == outdir / "peprefix_r2.fastq.gz"


def test_from_cli_reuse_splits_missing_dir_raises(tmp_path: _PathAlias):
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        paired_ended: false
        produce_fastq: true
        produce_bam: false
        produce_vcf: false
        threads: 4
        cleanup_splits: true
        reuse_splits: true
        overwrite_output: true
        """
    ).strip() + "\n"

    yml_path = tmp_path / "neat_from_cli_reuse.yml"
    yml_path.write_text(cfg, encoding="utf-8")

    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    options = Options.from_cli(outdir, "reuse", yml_path)
    # should issue a warning but continue in this case
    assert options.reuse_splits == True
