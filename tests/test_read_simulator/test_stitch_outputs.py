"""
Unit tests for neat/read_simulator/utils/stitch_outputs.py
"""
import gzip
import io
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch, call

import numpy as np
import pytest

from neat.read_simulator.utils.stitch_outputs import concat, merge_vcfs, merge_bam, main
from neat.variants import SingleNucleotideVariant
from neat.variants.contig_variants import ContigVariants


# Helpers

def _write_gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def _make_ofw(tmp_path: Path, vcf_path: Path = None):
    """
    Minimal OutputFileWriter stand-in backed by real StringIO / file handles.
    """
    ofw = SimpleNamespace()
    ofw.tmp_dir = tmp_path
    ofw.fq1 = tmp_path / "out.fq1.gz"
    ofw.fq2 = tmp_path / "out.fq2.gz"
    ofw.bam  = tmp_path / "out.bam"

    if vcf_path is None:
        vcf_path = tmp_path / "out.vcf.gz"
    ofw.vcf = vcf_path

    # Writable StringIO acts as the destination handle for concat/merge_vcfs
    ofw._fq1_buf = io.StringIO()
    ofw._fq2_buf = io.StringIO()
    ofw._vcf_buf = io.StringIO()

    ofw.files_to_write = {
        ofw.fq1: ofw._fq1_buf,
        ofw.fq2: ofw._fq2_buf,
        ofw.vcf: ofw._vcf_buf,
    }
    return ofw


# concat (byte-level gzip concatenation — concatenated gzip streams are valid gzip)

def test_concat_single_file(tmp_path):
    src = _write_gz(tmp_path / "a.gz", "hello\n")
    dest = tmp_path / "out.gz"
    concat([src], dest)
    with gzip.open(dest, "rt") as f:
        assert f.read() == "hello\n"


def test_concat_multiple_files(tmp_path):
    a = _write_gz(tmp_path / "a.gz", "line1\n")
    b = _write_gz(tmp_path / "b.gz", "line2\n")
    dest = tmp_path / "out.gz"
    concat([a, b], dest)
    with gzip.open(dest, "rt") as f:
        assert f.read() == "line1\nline2\n"


def test_concat_empty_list(tmp_path):
    """Empty input list: concat is a no-op and dest is not created."""
    dest = tmp_path / "out.gz"
    concat([], dest)
    assert not dest.exists()


def test_concat_preserves_content(tmp_path):
    content = "ACGT\nACGT\nACGT\n"
    src = _write_gz(tmp_path / "reads.gz", content)
    dest = tmp_path / "out.gz"
    concat([src], dest)
    with gzip.open(dest, "rt") as f:
        assert f.read() == content


def test_concat_order_is_preserved(tmp_path):
    files = [_write_gz(tmp_path / f"{i}.gz", f"chunk{i}\n") for i in range(5)]
    dest = tmp_path / "out.gz"
    concat(files, dest)
    with gzip.open(dest, "rt") as f:
        result = f.read()
    positions = [result.index(f"chunk{i}") for i in range(5)]
    assert positions == sorted(positions)


# merge_vcfs

def test_merge_vcfs_skips_comment_lines(tmp_path):
    vcf_text = "##header line\n#CHROM\tPOS\n1\t100\tA\tT\n"
    vcf = _write_gz(tmp_path / "v.vcf.gz", vcf_text)
    ofw = _make_ofw(tmp_path)
    merge_vcfs([vcf], ofw)
    result = ofw._vcf_buf.getvalue()
    assert "##header line" not in result
    assert "#CHROM" not in result
    assert "1\t100\tA\tT" in result


def test_merge_vcfs_multiple_files(tmp_path):
    v1 = _write_gz(tmp_path / "v1.vcf.gz", "##h\n1\t10\tA\tT\n")
    v2 = _write_gz(tmp_path / "v2.vcf.gz", "##h\n2\t20\tC\tG\n")
    ofw = _make_ofw(tmp_path)
    merge_vcfs([v1, v2], ofw)
    result = ofw._vcf_buf.getvalue()
    assert "1\t10\tA\tT" in result
    assert "2\t20\tC\tG" in result


def test_merge_vcfs_empty_list(tmp_path):
    ofw = _make_ofw(tmp_path)
    merge_vcfs([], ofw)
    assert ofw._vcf_buf.getvalue() == ""


def test_merge_vcfs_only_comments_produces_no_output(tmp_path):
    vcf = _write_gz(tmp_path / "v.vcf.gz", "##header\n#CHROM\n")
    ofw = _make_ofw(tmp_path)
    merge_vcfs([vcf], ofw)
    assert ofw._vcf_buf.getvalue() == ""


def test_merge_vcfs_preserves_data_line_order(tmp_path):
    lines = [f"chr1\t{i}\t.\tA\tT\t.\t.\t.\n" for i in range(1, 6)]
    vcf = _write_gz(tmp_path / "v.vcf.gz", "".join(lines))
    ofw = _make_ofw(tmp_path)
    merge_vcfs([vcf], ofw)
    result = ofw._vcf_buf.getvalue().splitlines()
    positions = [int(r.split("\t")[1]) for r in result]
    assert positions == list(range(1, 6))


def test_merge_vcfs_dedup_removes_identical_lines(tmp_path):
    """Identical lines from two thread VCFs are collapsed to one (Issue #256)."""
    line = "chr1\t100\t.\tA\tT\t42\tPASS\t.\tGT\t0|1\n"
    v1 = _write_gz(tmp_path / "t0.vcf.gz", line)
    v2 = _write_gz(tmp_path / "t1.vcf.gz", line)
    ofw = _make_ofw(tmp_path)
    merge_vcfs([v1, v2], ofw)
    result = [l for l in ofw._vcf_buf.getvalue().splitlines() if l.strip()]
    assert len(result) == 1


def test_merge_vcfs_distinct_lines_are_all_kept(tmp_path):
    """Distinct lines from two threads both appear in merged output."""
    v1 = _write_gz(tmp_path / "t0.vcf.gz", "chr1\t100\t.\tA\tT\t42\tPASS\t.\tGT\t0|1\n")
    v2 = _write_gz(tmp_path / "t1.vcf.gz", "chr1\t200\t.\tC\tG\t42\tPASS\t.\tGT\t0|1\n")
    ofw = _make_ofw(tmp_path)
    merge_vcfs([v1, v2], ofw)
    result = [l for l in ofw._vcf_buf.getvalue().splitlines() if l.strip()]
    assert len(result) == 2


def test_merge_vcfs_partial_overlap_deduped(tmp_path):
    """Three lines total, two of which are identical: result has two unique lines."""
    line_a = "chr1\t100\t.\tA\tT\t42\tPASS\t.\tGT\t0|1\n"
    line_b = "chr1\t200\t.\tC\tG\t42\tPASS\t.\tGT\t0|1\n"
    v1 = _write_gz(tmp_path / "t0.vcf.gz", line_a + line_b)
    v2 = _write_gz(tmp_path / "t1.vcf.gz", line_a)
    ofw = _make_ofw(tmp_path)
    merge_vcfs([v1, v2], ofw)
    result = [l for l in ofw._vcf_buf.getvalue().splitlines() if l.strip()]
    assert len(result) == 2


# find_dups (ContigVariants deduplication, Issue #256)

def test_find_dups_same_alt_different_genotype_rejected(tmp_path):
    """Same position + same ALT is a duplicate regardless of genotype."""
    cv = ContigVariants()
    v1 = SingleNucleotideVariant(10, "T", np.array([1, 0]), 40)
    v2 = SingleNucleotideVariant(10, "T", np.array([0, 1]), 40)
    cv.add_variant(v1)
    assert cv.add_variant(v2) == 1


def test_find_dups_different_alt_same_position_accepted(tmp_path):
    """Two SNVs at the same position with different ALTs are not duplicates."""
    cv = ContigVariants()
    v1 = SingleNucleotideVariant(10, "T", np.array([1, 0]), 40)
    v2 = SingleNucleotideVariant(10, "G", np.array([0, 1]), 40)
    cv.add_variant(v1)
    assert cv.add_variant(v2) == 0
    assert len(cv.contig_variants[10]) == 2


def test_find_dups_exact_duplicate_rejected(tmp_path):
    """Exact duplicates (same position, ALT, and genotype) are rejected."""
    cv = ContigVariants()
    v1 = SingleNucleotideVariant(10, "T", np.array([0, 1]), 40)
    v2 = SingleNucleotideVariant(10, "T", np.array([0, 1]), 40)
    cv.add_variant(v1)
    assert cv.add_variant(v2) == 1


# merge_bam (now a raw BGZF concatenation via pysam.cat — see stitch_outputs)

def test_merge_bam_calls_pysam_cat(tmp_path):
    ofw = _make_ofw(tmp_path)
    bam_files = [tmp_path / f"{i}.bam" for i in range(3)]

    with patch("neat.read_simulator.utils.stitch_outputs.pysam") as mock_pysam:
        merge_bam(bam_files, ofw, threads=4)

    # Per-worker BAMs are coord-sorted and cover non-overlapping ranges, so a single
    # pysam.cat is sufficient — no merge/sort needed.
    mock_pysam.cat.assert_called_once()
    mock_pysam.merge.assert_not_called()
    mock_pysam.sort.assert_not_called()


def test_merge_bam_cat_uses_output_bam_path(tmp_path):
    ofw = _make_ofw(tmp_path)
    bam_files = [tmp_path / "a.bam"]

    with patch("neat.read_simulator.utils.stitch_outputs.pysam") as mock_pysam:
        merge_bam(bam_files, ofw, threads=2)

    cat_args = mock_pysam.cat.call_args[0]
    assert "-o" in cat_args
    assert str(ofw.bam) in cat_args


def test_merge_bam_passes_all_inputs_to_cat(tmp_path):
    """Cat is given every per-worker BAM as a positional argument, in order."""
    ofw = _make_ofw(tmp_path)
    bam_files = [tmp_path / f"{i}.bam" for i in range(600)]

    with patch("neat.read_simulator.utils.stitch_outputs.pysam") as mock_pysam:
        merge_bam(bam_files, ofw, threads=1)

    cat_args = mock_pysam.cat.call_args[0]
    # cat call shape: ("-o", out_path, *input_paths)
    passed_inputs = [a for a in cat_args if a not in ("-o", str(ofw.bam))]
    assert passed_inputs == [str(p) for p in bam_files]


# main

def _file_dict(fq1=None, fq2=None, vcf=None, bam=None):
    return {"fq1": fq1, "fq2": fq2, "vcf": vcf, "bam": bam}


def test_main_fq1_only(tmp_path):
    ofw = _make_ofw(tmp_path)
    src = _write_gz(tmp_path / "chunk.fq1.gz", "@read1\nACGT\n+\nIIII\n")
    output_files = [(0, _file_dict(fq1=src))]
    main(ofw, output_files)
    # Byte-level concat writes raw bytes to ofw.fq1 (the path), not through the
    # StringIO handle. Read the actual file back.
    with gzip.open(ofw.fq1, "rt") as fh:
        assert "read1" in fh.read()


def test_main_fq1_and_fq2(tmp_path):
    ofw = _make_ofw(tmp_path)
    src1 = _write_gz(tmp_path / "c.fq1.gz", "@r1\nACGT\n+\nIIII\n")
    src2 = _write_gz(tmp_path / "c.fq2.gz", "@r2\nTTGG\n+\nIIII\n")
    output_files = [(0, _file_dict(fq1=src1, fq2=src2))]
    main(ofw, output_files)
    with gzip.open(ofw.fq1, "rt") as fh:
        assert "r1" in fh.read()
    with gzip.open(ofw.fq2, "rt") as fh:
        assert "r2" in fh.read()


def test_main_vcf(tmp_path):
    ofw = _make_ofw(tmp_path)
    src = _write_gz(tmp_path / "chunk.vcf.gz", "##header\nchr1\t100\t.\tA\tT\n")
    output_files = [(0, _file_dict(vcf=src))]
    main(ofw, output_files)
    result = ofw._vcf_buf.getvalue()
    assert "chr1\t100" in result
    assert "##header" not in result


def test_main_none_files_not_concatenated(tmp_path):
    """None entries in the file dict should be skipped without error."""
    ofw = _make_ofw(tmp_path)
    output_files = [(0, _file_dict())]  # all None
    main(ofw, output_files)  # should not raise
    # No inputs → byte-level concat is a no-op and the destination paths are not
    # created. The VCF buffer remains untouched (no merge_vcfs run).
    assert not ofw.fq1.exists()
    assert not ofw.fq2.exists()
    assert ofw._vcf_buf.getvalue() == ""


def test_main_multiple_threads(tmp_path):
    ofw = _make_ofw(tmp_path)
    chunks = [
        (i, _file_dict(fq1=_write_gz(tmp_path / f"c{i}.fq1.gz", f"chunk{i}\n")))
        for i in range(3)
    ]
    main(ofw, chunks)
    with gzip.open(ofw.fq1, "rt") as fh:
        result = fh.read()
    for i in range(3):
        assert f"chunk{i}" in result


def test_main_bam_calls_merge_bam(tmp_path):
    ofw = _make_ofw(tmp_path)
    bam_file = tmp_path / "chunk.bam"
    output_files = [(0, _file_dict(bam=bam_file))]

    with patch("neat.read_simulator.utils.stitch_outputs.merge_bam") as mock_merge:
        main(ofw, output_files, threads=2)

    mock_merge.assert_called_once_with([bam_file], ofw, 2)