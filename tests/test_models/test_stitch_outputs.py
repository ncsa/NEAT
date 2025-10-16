"""
Unit tests for helper functions in the stitch_outputs module
"""

from pathlib import Path

from Bio import bgzf

from neat.read_simulator.utils.stitch_outputs import concat


def test_concat_joins_files_in_order(tmp_path: Path) -> None:
    """Verify that concat writes the exact bytewise concatenation of inputs."""
    # Prepare small input files
    f1 = tmp_path / "a.bin"
    with bgzf.BgzfWriter(f1, 'w')as f1_in:
        f1_in.write("hello\n")
    f2 = tmp_path / "b.bin"
    with bgzf.BgzfWriter(f2, 'w') as f2_in:
        f2_in.write("world\n")
    f3 = tmp_path / "c.bin"
    with bgzf.BgzfWriter(f3, 'w') as f3_in:
        f3_in.write("!!!")

    dest = tmp_path / "out.bin"
    dest_write = bgzf.BgzfWriter(dest, 'w')
    concat([f1, f2, f3], dest_write)
    dest_write.close()
    assert dest.exists()
    with bgzf.BgzfReader(dest) as read_dest:
        text = ""
        for line in read_dest:
            text += line
    assert text == "hello\nworld\n!!!"


def test_concat_noop_on_empty_list(tmp_path: Path) -> None:
    dest = tmp_path / "out_empty.bin"
    concat([], dest)
    # Should not create the destination file when given no inputs
    assert not dest.exists()