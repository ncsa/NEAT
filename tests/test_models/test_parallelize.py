"""
Unit tests for the parallel_runner module
"""

import sys
from pathlib import Path
from neat.read_simulator import parallel_runner as pr


def test_worker_creates_dir_and_runs(tmp_path: Path) -> None:
    """Ensure that worker creates the working directory and executes the command."""
    workdir = tmp_path / "job1"
    # benign python command
    cmd = [sys.executable, "-c", "print('ok')"]
    returned = pr.worker((cmd, workdir))
    assert returned == workdir
    assert workdir.is_dir()


def test_worker_raises_on_failure(tmp_path: Path) -> None:
    """Worker should wrap failing subprocess in RuntimeError."""
    workdir = tmp_path / "job2"
    cmd = [sys.executable, "-c", "import sys; sys.exit(2)"]
    try:
        pr.worker((cmd, workdir))
    except RuntimeError as e:
        assert "NEAT job failed" in str(e)
        assert workdir.is_dir()
    else:
        raise AssertionError("worker did not raise RuntimeError on failing command")


def test_module_exports_and_constants() -> None:
    """Basic smoke checks for public API surface in the module."""
    assert callable(pr.main)
    assert isinstance(pr.EXTENSIONS, list)
    assert set(pr.EXTENSIONS) >= {"gz", "fastq", "bam", "vcf"}