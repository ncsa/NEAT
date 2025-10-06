"""
Unit tests for helper functions in the stitch_outputs module
"""

from pathlib import Path
from neat.read_simulator.utils.stitch_outputs import concat


def test_concat_joins_files_in_order(tmp_path: Path) -> None:
    """Verify that concat writes the exact bytewise concatenation of inputs."""
    # Prepare small input files
    f1 = tmp_path / "a.bin"
    f1.write_bytes(b"hello\n")
    f2 = tmp_path / "b.bin"
    f2.write_bytes(b"world\n")
    f3 = tmp_path / "c.bin"
    f3.write_bytes(b"!!!")

    dest = tmp_path / "out.bin"
    concat([f1, f2, f3], dest)

    assert dest.exists()
    assert dest.read_bytes() == b"hello\nworld\n!!!"


def test_concat_noop_on_empty_list(tmp_path: Path) -> None:
    dest = tmp_path / "out_empty.bin"
    concat([], dest)
    # Should not create the destination file when given no inputs
    assert not dest.exists()