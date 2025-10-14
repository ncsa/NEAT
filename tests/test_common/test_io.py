import os
import stat
import gzip
from pathlib import Path

import pytest
from Bio import bgzf

from neat.common.io import is_compressed, open_input, open_output, validate_input_path, validate_output_path


def test_is_compressed_plain_and_gz(tmp_path: Path):
    plain = tmp_path / "file.txt"
    plain.write_text("hello", encoding="utf-8")
    gz = tmp_path / "file.txt.gz"
    with bgzf.BgzfWriter(gz, "wt") as f:
        f.write("hello")

    assert is_compressed(plain) is False
    assert is_compressed(gz) is True


def test_open_input_reads_plain_and_gz(tmp_path: Path):
    # Plain
    plain = tmp_path / "a.txt"
    plain.write_text("abc\n", encoding="utf-8")
    with open_input(plain) as fh:
        assert fh.read() == "abc\n"

    # Gz
    gz = tmp_path / "b.txt.gz"
    with bgzf.BgzfWriter(gz, "wt") as fh:
        fh.write("xyz\n")
    out = ""
    with open_input(gz) as fh:
        for line in fh:
            out += line
    assert out == "xyz\n"


def test_open_output_creates_dirs_and_writes_plain(tmp_path: Path):
    out = tmp_path / "nested/dir/out.txt"
    with open_output(out) as fh:
        fh.write("data")
    assert out.exists()
    assert out.read_text(encoding="utf-8") == "data"


def test_open_output_handles_gz_and_bgz(tmp_path: Path):
    # .gz
    out_gz = tmp_path / "out1.fq.gz"
    with open_output(out_gz) as fh:
        fh.write("line1\nline2\n")
    # Read back with gzip
    with gzip.open(out_gz, "rt", encoding="utf-8") as fh:
        assert fh.read() == "line1\nline2\n"

    # .bgz
    out_bgz = tmp_path / "out2.fq.bgz"
    with open_output(out_bgz) as fh:
        fh.write("a\nb\n")
    with gzip.open(out_bgz, "rt", encoding="utf-8") as fh:
        assert fh.read() == "a\nb\n"


def test_open_output_xt_mode_existing_bgz_exits(tmp_path: Path):
    p = tmp_path / "exists.bgz"
    # create file first
    with open_output(p) as fh:
        fh.write("x")
    # now opening with xt should exit with code 3
    with pytest.raises(SystemExit) as ei:
        with open_output(p, mode="xt"):
            pass
    assert ei.value.code == 3


def test_validate_input_path_errors(tmp_path: Path):
    # non-existent
    with pytest.raises(SystemExit) as ei:
        validate_input_path(tmp_path / "missing.txt")
    assert ei.value.code == 5

    # empty file
    empty = tmp_path / "empty.txt"
    empty.touch()
    with pytest.raises(SystemExit) as ei:
        validate_input_path(empty)
    assert ei.value.code == 7

    # no read permission
    f = tmp_path / "noread.txt"
    f.write_text("x", encoding="utf-8")
    # remove read permission for owner
    current = f.stat().st_mode
    f.chmod(current & ~stat.S_IRUSR & ~stat.S_IRGRP & ~stat.S_IROTH)
    try:
        with pytest.raises(SystemExit) as ei:
            validate_input_path(f)
        assert ei.value.code == 9
    finally:
        # restore permission so tmp cleanup works
        f.chmod(current)


def test_validate_output_path_file_and_dir(tmp_path: Path):
    # file exists, no overwrite
    f = tmp_path / "out.txt"
    f.write_text("x", encoding="utf-8")
    with pytest.raises(SystemExit) as ei:
        validate_output_path(f, is_file=True, overwrite=False)
    assert ei.value.code == 3

    # overwrite allowed should not raise
    validate_output_path(f, is_file=True, overwrite=True)

    # directory exists but not writable
    d = tmp_path / "dir"
    d.mkdir()
    # remove write permission
    mode = d.stat().st_mode
    d.chmod(mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH)
    try:
        with pytest.raises(SystemExit) as ei:
            validate_output_path(d, is_file=False)
        assert ei.value.code == 11
    finally:
        d.chmod(mode)

    # nonexistent directory should be created
    d2 = tmp_path / "new/dir"
    validate_output_path(d2, is_file=False)
    assert d2.exists() and d2.is_dir()
