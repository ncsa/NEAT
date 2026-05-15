"""
Unit tests for the split_inputs module of the parallel read simulator
"""
import gzip

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.read_simulator.utils.options import Options
from neat.read_simulator.utils.split_inputs import (
    chunk_record,
    disk_bytes_free,
    main,
    print_stderr,
    write_fasta,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_options(tmp_path, mode="contig", read_len=50, block_size=200):
    opts = Options(rng_seed=0)
    opts.read_len = read_len
    opts.parallel_mode = mode
    opts.parallel_block_size = block_size
    opts.splits_dir = tmp_path / "splits"
    opts.splits_dir.mkdir(exist_ok=True)
    opts.reference = tmp_path / "ref.fa"
    return opts


def _make_ref_file(path, *contigs: tuple[str, str]):
    with open(path, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n{seq}\n")


def _read_fasta_gz(path) -> tuple[str, str]:
    """Return (header_line, sequence) from a bgzf/gzip FASTA."""
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    header = lines[0]
    seq = "".join(lines[1:])
    return header, seq


# ===========================================================================
# chunk_record
# ===========================================================================

def test_chunk_record_overlaps() -> None:
    """Existing behaviour: overlapping chunks of a 12-base sequence."""
    rec = Seq("ACGTACGTACGT")
    chunks = list(chunk_record(rec, 5, 2))
    assert len(chunks) == 4
    lengths = [len(seq) for _, seq in chunks]
    assert lengths == [5, 5, 5, 3]
    ids = [cid for cid, _ in chunks]
    seqs = [seq for _, seq in chunks]
    assert ids == [0, 3, 6, 9]
    assert seqs == [Seq("ACGTA"), Seq("TACGT"), Seq("GTACG"), Seq("CGT")]


def test_chunk_record_sequence_shorter_than_chunk():
    """Sequence shorter than chunk_len produces one chunk covering the whole sequence."""
    rec = Seq("ACGT")
    chunks = list(chunk_record(rec, 100, 10))
    assert len(chunks) == 1
    start, seq = chunks[0]
    assert start == 0
    assert seq == rec


def test_chunk_record_sequence_exact_chunk_length():
    """Sequence exactly equal to chunk_len produces exactly one chunk."""
    rec = Seq("ACGTACGT")
    chunks = list(chunk_record(rec, 8, 2))
    assert len(chunks) == 1
    assert chunks[0] == (0, rec)


def test_chunk_record_no_overlap():
    """With overlap=0 chunks are non-overlapping and contiguous."""
    rec = Seq("ACGTACGT")
    chunks = list(chunk_record(rec, 4, 0))
    assert len(chunks) == 2
    assert chunks[0] == (0, Seq("ACGT"))
    assert chunks[1] == (4, Seq("ACGT"))


def test_chunk_record_start_positions_are_correct():
    """Start of each chunk equals end_of_previous - overlap."""
    rec = Seq("A" * 100)
    chunk_len, overlap = 30, 10
    chunks = list(chunk_record(rec, chunk_len, overlap))
    for i in range(1, len(chunks)):
        prev_end = chunks[i - 1][0] + len(chunks[i - 1][1])
        assert chunks[i][0] == prev_end - overlap


def test_chunk_record_covers_full_sequence():
    """The union of all chunks covers every position in the sequence."""
    seq = "ACGT" * 25  # 100 bp
    rec = Seq(seq)
    chunks = list(chunk_record(rec, 30, 5))
    covered = set()
    for start, subseq in chunks:
        for i in range(len(subseq)):
            covered.add(start + i)
    assert covered == set(range(len(seq)))


# ===========================================================================
# print_stderr
# ===========================================================================

def test_print_stderr_no_exit():
    """Calling with exit_=False should not raise."""
    print_stderr("test message", exit_=False)  # no exception


def test_print_stderr_exit():
    """Calling with exit_=True should raise SystemExit."""
    with pytest.raises(SystemExit):
        print_stderr("fatal error", exit_=True)


# ===========================================================================
# disk_bytes_free
# ===========================================================================

def test_disk_bytes_free_returns_positive(tmp_path):
    result = disk_bytes_free(tmp_path)
    assert isinstance(result, int)
    assert result > 0


# ===========================================================================
# write_fasta
# ===========================================================================

def test_write_fasta_creates_file(tmp_path):
    out = tmp_path / "test.fa.gz"
    write_fasta("chr1", Seq("ACGTACGT"), out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_write_fasta_header(tmp_path):
    out = tmp_path / "test.fa.gz"
    write_fasta("mycontig", Seq("ACGT"), out)
    header, _ = _read_fasta_gz(out)
    assert header == ">mycontig"


def test_write_fasta_sequence_content(tmp_path):
    seq = "ACGT" * 10
    out = tmp_path / "test.fa.gz"
    write_fasta("chr1", Seq(seq), out)
    _, recovered = _read_fasta_gz(out)
    assert recovered == seq


def test_write_fasta_line_wrapping_default(tmp_path):
    """Default width=80 wraps sequences longer than 80 bases."""
    seq = "A" * 200
    out = tmp_path / "test.fa.gz"
    write_fasta("chr1", Seq(seq), out)
    with gzip.open(out, "rt") as fh:
        lines = fh.read().splitlines()
    seq_lines = [l for l in lines if not l.startswith(">")]
    assert all(len(l) <= 80 for l in seq_lines)
    assert "".join(seq_lines) == seq


def test_write_fasta_custom_width(tmp_path):
    seq = "A" * 50
    out = tmp_path / "test.fa.gz"
    write_fasta("chr1", Seq(seq), out, width=10)
    with gzip.open(out, "rt") as fh:
        lines = fh.read().splitlines()
    seq_lines = [l for l in lines if not l.startswith(">")]
    assert all(len(l) <= 10 for l in seq_lines)


# ===========================================================================
# main — contig mode
# ===========================================================================

def test_main_contig_mode_one_file_per_contig(tmp_path):
    opts = _make_options(tmp_path, mode="contig")
    _make_ref_file(opts.reference, ("chr1", "ACGT" * 50), ("chr2", "TTGG" * 30))
    result, written, ref_lens = main(opts)
    assert written == 2
    assert set(result.keys()) == {"chr1", "chr2"}
    assert ref_lens == {"chr1": 200, "chr2": 120}


def test_main_contig_mode_dict_keyed_by_full_span(tmp_path):
    opts = _make_options(tmp_path, mode="contig")
    _make_ref_file(opts.reference, ("chr1", "ACGT" * 25))  # 100 bp
    result, _, _ = main(opts)
    keys = list(result["chr1"].keys())
    assert len(keys) == 1
    assert keys[0] == (0, 100)


def test_main_contig_mode_files_exist(tmp_path):
    opts = _make_options(tmp_path, mode="contig")
    _make_ref_file(opts.reference, ("chr1", "ACGT" * 10), ("chr2", "TTGG" * 10))
    result, _, _ = main(opts)
    for contig, spans in result.items():
        for path in spans.values():
            assert path.exists(), f"Expected {path} to exist"


def test_main_contig_mode_filename_contains_contig(tmp_path):
    opts = _make_options(tmp_path, mode="contig")
    _make_ref_file(opts.reference, ("mychr", "ACGT" * 10))
    result, _, _ = main(opts)
    path = list(result["mychr"].values())[0]
    assert "mychr" in path.name


def test_main_contig_mode_filename_zero_padded(tmp_path):
    opts = _make_options(tmp_path, mode="contig")
    _make_ref_file(opts.reference, ("chr1", "A" * 40), ("chr2", "T" * 40))
    result, _, _ = main(opts)
    names = sorted(p.name for spans in result.values() for p in spans.values())
    # First file index should be zero-padded to width 10
    assert names[0].startswith("0000000001")
    assert names[1].startswith("0000000002")


def test_main_contig_mode_fasta_content(tmp_path):
    """Files written in contig mode should contain the correct uppercased sequence."""
    seq = "acgtACGT" * 5  # mixed case
    opts = _make_options(tmp_path, mode="contig")
    _make_ref_file(opts.reference, ("chr1", seq))
    result, _, _ = main(opts)
    path = list(result["chr1"].values())[0]
    _, recovered = _read_fasta_gz(path)
    assert recovered == seq.upper()


# ===========================================================================
# main — block (size) mode
# ===========================================================================

def test_main_block_mode_produces_multiple_chunks(tmp_path):
    """A sequence much longer than block_size should produce multiple chunks."""
    opts = _make_options(tmp_path, mode="size", read_len=10, block_size=50)
    seq = "ACGT" * 50  # 200 bp — should produce several 50 bp blocks
    _make_ref_file(opts.reference, ("chr1", seq))
    result, written, _ = main(opts)
    assert written > 1
    assert len(result["chr1"]) > 1


def test_main_block_mode_written_count_matches_dict(tmp_path):
    opts = _make_options(tmp_path, mode="size", read_len=10, block_size=50)
    _make_ref_file(opts.reference, ("chr1", "A" * 200), ("chr2", "T" * 200))
    result, written, _ = main(opts)
    total_spans = sum(len(spans) for spans in result.values())
    assert written == total_spans


def test_main_block_mode_files_exist(tmp_path):
    opts = _make_options(tmp_path, mode="size", read_len=10, block_size=50)
    _make_ref_file(opts.reference, ("chr1", "ACGT" * 50))
    result, _, _ = main(opts)
    for spans in result.values():
        for path in spans.values():
            assert path.exists()


def test_main_block_mode_span_keys_cover_sequence(tmp_path):
    """Span keys should cover positions 0 through len(seq), with overlaps."""
    opts = _make_options(tmp_path, mode="size", read_len=10, block_size=50)
    seq = "A" * 200
    _make_ref_file(opts.reference, ("chr1", seq))
    result, _, _ = main(opts)
    spans = sorted(result["chr1"].keys())
    assert spans[0][0] == 0
    # Last span should end at or cover the full sequence length
    assert spans[-1][1] >= len(seq)


def test_main_block_mode_sequence_shorter_than_block(tmp_path):
    """A short sequence produces exactly one chunk even in block mode."""
    opts = _make_options(tmp_path, mode="size", read_len=10, block_size=500)
    _make_ref_file(opts.reference, ("chr1", "ACGT" * 10))  # 40 bp < 500
    result, written, _ = main(opts)
    assert written == 1
    assert len(result["chr1"]) == 1