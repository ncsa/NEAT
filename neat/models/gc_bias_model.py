"""
GC-bias model for NEAT.
"""

import logging
import pickle
import numpy as np
from typing import Union
from pathlib import Path

_LOG = logging.getLogger(__name__)

class GCBiasModel:
    """
    Per-GC-content weight table used to bias fragment retention during read simulation.
    Weights are stored as 101 values indexed by integer GC percentage (0-100%).
    """

    def __init__(self, weights: Union[list, np.ndarray], window_size: int):
        if len(weights) != 101:
            raise ValueError("GC bias weights must have exactly 101 elements.")
        if window_size <= 0:
            raise ValueError("GC bias window size must be positive.")
        
        self.weights = np.array(weights, dtype=float)
        if np.any(self.weights < 0):
            raise ValueError("GC bias weights must be non-negative.")
        if np.all(self.weights == 0):
            raise ValueError("GC bias model must contain at least one positive weight.")
            
        self.window_size = window_size
        self.is_uniform = np.all(self.weights == self.weights[0])
        self.max_weight = np.max(self.weights)

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> 'GCBiasModel':
        """
        Loads GC bias model from a pickle file.
        The pickle file is expected to contain [GC_SCALE_COUNT, GC_SCALE_VAL]
        where GC_SCALE_COUNT[-1] is the window size and GC_SCALE_VAL is the weights.
        """
        with open(path, 'rb') as f:
            data = pickle.load(f)
            
        if isinstance(data, list) and len(data) == 2:
            gc_scale_count, weights = data
            window_size = gc_scale_count[-1]
            return cls(weights, window_size)
        else:
            raise ValueError(f"Unexpected data format in GC bias model file: {path}")

    def save(self, path: Union[str, Path]):
        """
        Saves GC bias model to a pickle file in the format compatible with NEAT 2.1.
        """
        gc_scale_count = list(range(1, 101)) + [self.window_size]
        data = [gc_scale_count, self.weights.tolist()]
        with open(path, 'wb') as f:
            pickle.dump(data, f)

    def get_weight(self, gc_fraction: float) -> float:
        """
        Returns the weight for a given GC fraction.
        """
        if self.is_uniform:
            return self.weights[0]
            
        index = int(round(gc_fraction * 100))
        index = max(0, min(100, index))
        return self.weights[index]

    def get_weight_for_sequence(self, sequence: str) -> float:
        """
        Calculates GC fraction of the sequence and returns the corresponding weight.
        """
        if self.is_uniform:
            return self.weights[0]
            
        called_bases = 0
        gc_count = 0
        for base in sequence.upper():
            if base in 'GC':
                gc_count += 1
                called_bases += 1
            elif base in 'AT':
                called_bases += 1
        
        if called_bases == 0:
            return 1.0 # Neutral weight for all-N sequences
            
        gc_fraction = gc_count / called_bases
        return self.get_weight(gc_fraction)

def get_uniform_gc_model(window_size: int = 100) -> GCBiasModel:
    """
    Returns a uniform GC bias model (no bias).
    """
    return GCBiasModel([1.0] * 101, window_size)
