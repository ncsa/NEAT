"""
Unit tests for the split_inputs module of the parallel read simulator
"""

from neat.parallel_read_simulator.split_inputs import SimpleRecord, chunk_record


def test_chunk_record_overlaps() -> None:
    """Ensure that chunk_record yields overlapping chunks with the correct ids."""
    # A simple sequence of 12 bases to chunk into length 5 with overlap 2
    rec = SimpleRecord("contig", "ACGTACGTACGT")
    chunks = list(chunk_record(rec, 5, 2))
    # Should produce four chunks: positions [0:5], [3:8], [6:11], [9:12]
    assert len(chunks) == 4
    # Check that chunk ids are sequential and lengths match expected slices
    lengths = [len(r.seq) for r, _ in chunks]
    assert lengths == [5, 5, 5, 3]
    ids = [cid for _, cid in chunks]
    assert ids == [1, 2, 3, 4]