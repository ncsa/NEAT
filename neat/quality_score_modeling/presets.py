"""
Named quality-bin presets for common Illumina instruments.

Each preset maps to the discrete Phred score levels that the instrument emits.
Pass the preset name to ``neat model-qual-score --quality-preset`` instead of
specifying raw bin values with ``-Q``.

References:
  NovaSeq 6000 / NovaSeq X: Illumina RTA3 4-level binning (Q2, Q12, Q23, Q37)
  NextSeq 2000: 4-level binning matching NovaSeq X (Q2, Q12, Q26, Q37)
  NextSeq 500 / MiniSeq: 5-level binning (Q2, Q12, Q23, Q27, Q37)
"""

__all__ = ["QUALITY_PRESETS"]

QUALITY_PRESETS: dict[str, list[int]] = {
    "novaseq":     [2, 12, 23, 37],
    "nextseq2000": [2, 12, 26, 37],
    "nextseq500":  [2, 12, 23, 27, 37],
}