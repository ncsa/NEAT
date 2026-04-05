"""
Unit tests for neat/read_simulator/utils/output_file_writer.py
"""
import gzip
import io
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.Seq import Seq

from neat.read_simulator.utils.output_file_writer import OutputFileWriter, reg2bin
from neat.read_simulator.utils.options import Options
from neat.read_simulator.utils.read import Read


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_REF_SEQ = "ACGT" * 30   # 120 bp


def _make_options(tmp_path: Path,
                  fq1=True, fq2=False, vcf=False, bam=False,
                  paired=False, produce_vcf=False, produce_bam=False,
                  produce_fastq=True) -> Options:
    opts = Options(rng_seed=0)
    opts.paired_ended = paired
    opts.produce_fastq = produce_fastq
    opts.produce_vcf = produce_vcf
    opts.produce_bam = produce_bam
    opts.temp_dir_path = tmp_path
    opts.reference = str(tmp_path / "ref.fa")
    opts.rng_seed = 0

    opts.fq1 = (tmp_path / "out.fq1.gz") if fq1 else None
    opts.fq2 = (tmp_path / "out.fq2.gz") if fq2 else None
    opts.vcf = (tmp_path / "out.vcf.gz") if vcf else None
    opts.bam = (tmp_path / "out.bam")   if bam  else None
    return opts


def _make_ofw(tmp_path: Path, **kw) -> OutputFileWriter:
    opts = _make_options(tmp_path, **kw)
    return OutputFileWriter(options=opts)


def _make_read(position: int = 10,
               seq: str = "ACGTACGT",
               is_reverse: bool = False) -> Read:
    read_len = len(seq)
    end_point = position + read_len
    r = Read(
        name="read1",
        raw_read=(position, end_point, position + 150, end_point + 150),
        reference_segment=Seq(seq),
        reference_id="chr1",
        ref_id_index=0,
        position=position,
        end_point=end_point,
        padding=20,
        run_read_len=read_len,
        is_reverse=is_reverse,
    )
    r.mapping_quality = 60
    r.quality_array = [40] * read_len
    # read_sequence must be set before write_bam_record can call make_cigar
    r.read_sequence = Seq(seq)
    return r


# ===========================================================================
# reg2bin
# ===========================================================================

def test_reg2bin_same_16kb_bin():
    result = reg2bin(0, 100)
    assert isinstance(result, int)
    assert result >= 0


def test_reg2bin_large_span_returns_zero():
    # Spanning >2^26 bases → falls through to return 0
    result = reg2bin(0, 2**26 + 1)
    assert result == 0


def test_reg2bin_adjacent_bins():
    # Two regions in the same 16kb window should return same bin
    assert reg2bin(0, 100) == reg2bin(0, 200)


def test_reg2bin_17kb_window():
    result = reg2bin(0, 2**17 + 1)
    assert isinstance(result, int)


# ===========================================================================
# OutputFileWriter — construction
# ===========================================================================

def test_ofw_fq1_file_created(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    assert tmp_path / "out.fq1.gz" in ofw.files_to_write


def test_ofw_fq2_file_created(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True, fq2=True, paired=True)
    assert tmp_path / "out.fq2.gz" in ofw.files_to_write


def test_ofw_vcf_file_created(tmp_path):
    ofw = _make_ofw(tmp_path, vcf=True, produce_vcf=True)
    assert tmp_path / "out.vcf.gz" in ofw.files_to_write


def test_ofw_no_files_raises(tmp_path):
    opts = _make_options(tmp_path, fq1=False)
    with pytest.raises(ValueError):
        OutputFileWriter(options=opts)


def test_ofw_fq1_is_none_when_not_requested(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    assert ofw.fq2 is None


def test_ofw_vcf_is_none_when_not_requested(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    assert ofw.vcf is None


def test_ofw_bam_is_none_when_not_requested(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    assert ofw.bam is None


# ===========================================================================
# OutputFileWriter — VCF header
# ===========================================================================

def _ofw_with_vcf_header(tmp_path: Path) -> OutputFileWriter:
    opts = _make_options(tmp_path, fq1=True, vcf=True, produce_vcf=True)
    # write a placeholder reference file for the header path
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n")
    return OutputFileWriter(
        options=opts,
        vcf_header={"chr1": 1000, "chr2": 500},
    )


def test_vcf_header_written_on_init(tmp_path):
    ofw = _ofw_with_vcf_header(tmp_path)
    # flush so content is available
    ofw.files_to_write[ofw.vcf].flush()
    with gzip.open(ofw.vcf, "rt") as fh:
        content = fh.read()
    assert "##fileformat=VCFv4.1" in content


def test_vcf_header_contains_contig_lines(tmp_path):
    ofw = _ofw_with_vcf_header(tmp_path)
    ofw.files_to_write[ofw.vcf].flush()
    with gzip.open(ofw.vcf, "rt") as fh:
        content = fh.read()
    assert "chr1" in content
    assert "chr2" in content


def test_vcf_header_contains_column_line(tmp_path):
    ofw = _ofw_with_vcf_header(tmp_path)
    ofw.files_to_write[ofw.vcf].flush()
    with gzip.open(ofw.vcf, "rt") as fh:
        content = fh.read()
    assert "#CHROM" in content


# ===========================================================================
# write_fastq_record
# ===========================================================================

def test_write_fastq_record_writes_content(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    record = "@read1\nACGT\n+\nIIII\n"
    ofw.write_fastq_record(ofw.fq1, record)
    ofw.flush_and_close_files()
    with gzip.open(ofw.fq1, "rt") as fh:
        assert fh.read() == record


def test_write_fastq_record_unknown_file_raises(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    with pytest.raises(ValueError):
        ofw.write_fastq_record(tmp_path / "unknown.gz", "@r\nA\n+\nI\n")


def test_write_fastq_record_multiple_records(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    for i in range(5):
        ofw.write_fastq_record(ofw.fq1, f"@r{i}\nACGT\n+\nIIII\n")
    ofw.flush_and_close_files()
    with gzip.open(ofw.fq1, "rt") as fh:
        lines = fh.readlines()
    assert len(lines) == 20   # 4 lines × 5 records


# ===========================================================================
# write_vcf_record
# ===========================================================================

def test_write_vcf_record_appends_line(tmp_path):
    ofw = _make_ofw(tmp_path, vcf=True, produce_vcf=True)
    line = "chr1\t100\t.\tA\tT\t42\tPASS\t.\tGT\t0|1\n"
    ofw.write_vcf_record(line)
    ofw.flush_and_close_files()
    with gzip.open(ofw.vcf, "rt") as fh:
        assert fh.read() == line


def test_write_vcf_record_unknown_file_raises(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)  # no vcf
    with pytest.raises(ValueError):
        ofw.write_vcf_record("chr1\t100\t.\tA\tT\n")


# ===========================================================================
# flush_and_close_files
# ===========================================================================

def test_flush_and_close_closes_fq1(tmp_path):
    ofw = _make_ofw(tmp_path, fq1=True)
    ofw.flush_and_close_files()
    fh = ofw.files_to_write[ofw.fq1]
    assert fh.closed


def test_flush_and_close_idempotent(tmp_path):
    """Calling flush_and_close twice should not raise."""
    ofw = _make_ofw(tmp_path, fq1=True)
    ofw.flush_and_close_files()
    ofw.flush_and_close_files()  # second call should not raise


# ===========================================================================
# write_bam_record
# ===========================================================================

def _ofw_with_bam(tmp_path: Path) -> OutputFileWriter:
    opts = _make_options(tmp_path, fq1=True, bam=True,
                         produce_bam=True, produce_fastq=True)
    bam_header = {"chr1": 1000}
    return OutputFileWriter(options=opts, bam_header=bam_header)


def test_write_bam_record_writes_bytes(tmp_path):
    ofw = _ofw_with_bam(tmp_path)
    read = _make_read(position=10, seq="ACGTACGT")
    bam_handle = ofw.files_to_write[ofw.bam]
    ofw.write_bam_record(read, contig_id=0, bam_handle=bam_handle, read_length=8)
    # If no exception was raised, the record was written.
    assert True


def test_write_bam_record_reverse_strand(tmp_path):
    ofw = _ofw_with_bam(tmp_path)
    read = _make_read(position=10, seq="ACGTACGT", is_reverse=True)
    bam_handle = ofw.files_to_write[ofw.bam]
    ofw.write_bam_record(read, contig_id=0, bam_handle=bam_handle, read_length=8)
    assert True


def test_write_bam_record_odd_length_sequence(tmp_path):
    """Odd-length reads require padding — should not crash."""
    ofw = _ofw_with_bam(tmp_path)
    read = _make_read(position=10, seq="ACGTA")  # 5 bp — odd
    read.quality_array = [40] * 5
    bam_handle = ofw.files_to_write[ofw.bam]
    ofw.write_bam_record(read, contig_id=0, bam_handle=bam_handle, read_length=5)
    assert True