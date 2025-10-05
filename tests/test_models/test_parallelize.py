"""
Unit tests for the argument parsing in the parallelize module
"""

from pathlib import Path
from neat.read_simulator.parallel_runner import parse_args


def test_parse_args_defaults() -> None:
    """Ensure that parse_args sets sensible defaults for optional flags."""
    args = parse_args(["config.yml"])
    # Positional argument is preserved as a Path
    assert args.config == Path("config.yml")
    # Optional arguments default to None prior to merging with config
    assert args.outdir is None
    assert args.by is None
    assert args.size is None
    assert args.cleanup_splits is None
    assert args.reuse_splits is None
    assert args.jobs is None
    assert args.neat_cmd is None
    assert args.samtools is None
    assert args.final_prefix is None