from neat.read_simulator.utils.options import Options

from pathlib import Path as _PathAlias
import logging as _logging
import numpy as _np
import textwrap as _textwrap
import pytest as _pytest


def _project_root() -> _PathAlias:
    return _PathAlias(__file__).resolve().parents[2]


@_pytest.fixture(autouse=True)
def _isolate_neat_logging():
    """
    Prevent flaky 'ValueError: I/O operation on closed file' logging errors under pytest.
    """
    # Clear handlers on NEAT and all child loggers
    for name, logger in list(_logging.Logger.manager.loggerDict.items()):
        if name == "neat" or name.startswith("neat."):
            if isinstance(logger, _logging.Logger):
                for h in list(logger.handlers):
                    logger.removeHandler(h)
                    try:
                        h.close()
                    except Exception:
                        pass
                logger.handlers.clear()
                logger.propagate = True # child loggers will propagate to 'neat'

    neat_logger = _logging.getLogger("neat")
    neat_logger.handlers.clear()
    neat_logger.addHandler(_logging.NullHandler())
    neat_logger.propagate = False # stop at 'neat' (do not reach root)

    yield

    # Rremove NullHandler
    for h in list(neat_logger.handlers):
        neat_logger.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass
    neat_logger.handlers.clear()
    neat_logger.propagate = True


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

        parallel_mode: size
        parallel_block_size: 500000
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

    assert opts.reference == _project_root() / "data" / "H1N1.fa"
    assert opts.read_len == 75
    assert opts.coverage == 5
    assert opts.ploidy == 2
    assert opts.rng_seed == 42

    assert opts.output_dir == outdir
    assert opts.output_prefix == "fromcli"
    assert opts.fq1 == outdir / "fromcli.fastq.gz"
    assert opts.fq2 is None
    assert opts.bam is None
    assert opts.vcf is None

    assert opts.threads == 2
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

        parallel_mode: contig
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
        parallel_mode: size
        parallel_block_size: 500000
        cleanup_splits: true
        reuse_splits: true
        overwrite_output: true
        """
    ).strip() + "\n"

    yml_path = tmp_path / "neat_from_cli_reuse.yml"
    yml_path.write_text(cfg, encoding="utf-8")

    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    with _pytest.raises(FileNotFoundError, match=r"reuse_splits=True"):
        Options.from_cli(outdir, "reuse", yml_path)