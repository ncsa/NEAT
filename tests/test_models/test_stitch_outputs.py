"""
Unit tests for helper functions in the stitch_outputs module
"""

from pathlib import Path
from neat.parallel_read_simulator.stitch_outputs import natural_key


def test_natural_key_sorting() -> None:
    """Verify that natural_key sorts numeric suffixes in human order."""
    names = ["sample10.txt", "sample2.txt", "sample1.txt", "sample1a.txt"]
    paths = [Path(n) for n in names]
    sorted_paths = sorted(paths, key=natural_key)
    assert [p.name for p in sorted_paths] == ["sample1.txt", "sample1a.txt", "sample2.txt", "sample10.txt"]