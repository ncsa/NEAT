"""
Unit tests for the split_inputs module of the parallel read simulator
"""
from Bio.Seq import Seq

from neat.read_simulator.utils.split_inputs import chunk_record

def test_chunk_record_overlaps() -> None:
    """Ensure that chunk_record yields overlapping chunks with the correct ids."""
    # A simple sequence of 12 bases to chunk into length 5 with overlap 2
    rec = Seq("ACGTACGTACGT")
    chunks = list(chunk_record(rec, 5, 2))
    # Should produce four chunks: positions [0:5], [3:8], [6:11], [9:12]
    assert len(chunks) == 4
    # Check that chunk ids are sequential and lengths match expected slices
    lengths = [len(seq) for _, seq in chunks]
    assert lengths == [5, 5, 5, 3]
    ids = [cid for cid, _ in chunks]
    seqs = [seq for _, seq in chunks]
    assert ids == [0, 3, 6, 9]
    assert seqs == [Seq('ACGTA'), Seq('TACGT'), Seq('GTACG'), Seq('CGT')]