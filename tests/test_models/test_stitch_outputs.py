"""
Unit tests for helper functions in the stitch_outputs module
"""

from pathlib import Path

from Bio import bgzf

from neat.read_simulator.utils.stitch_outputs import concat


    # test_concat_joins_files_in_order removed:
    # covered by test_read_simulator/test_stitch_outputs.py::test_concat_multiple_files
    # and test_concat_preserves_content_order


def test_concat_noop_on_empty_list(tmp_path: Path) -> None:
    dest = tmp_path / "out_empty.bin"
    concat([], dest)
    # Should not create the destination file when given no inputs
    assert not dest.exists()