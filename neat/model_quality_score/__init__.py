"""
Submodule to build quality score models for NEAT that wraps NEAT’s
existing sequencing error model and traditional quality model, with the
option to construct a Markov chain–based quality model instead
"""

__all__ = ["model_qual_score_runner"]

from .runner import model_qual_score_runner